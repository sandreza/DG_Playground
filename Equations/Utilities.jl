"""
    Experimenting with macros for defining core types.
    This code is based on the UFL Julia macros
    written by Andrew Monteith:
    https://github.com/AndrewMonteith/UFL/blob/master/src/uflmacros.jl
"""

function get_opt_type(type)
    types = if isa(type, Symbol)
        [type]
    else
        union_types = type.args

        union_types[2:end]
    end
    return :(Union{$(types...), Nothing})
end

macro opt_t(type)
    return esc(get_opt_type(type))
end

macro opt(e)
    # Takes a variable declaration and x::T and transforms it to x::Union{T, Nothing}=nothing
    # If e is a x::Union{T1, ..., TN} it becomes x::Union{T1, ..., TN, Nothing}=nothing
    sym, types = e.args[1], get_opt_type(e.args[2])
    ex = :($sym::$types=nothing)
    ex.head = :kw
    return esc(ex)
end

export core_type, attach_hash_operators

field(sym::Symbol, t) = Expr(:(::), sym, t)

fields = Dict(
    :shape => (type=:DimensionTuple, default_val=()),
    :operands => (type=VarTuple{AbstractExpr}, default_val=()),
)
        
function find_field(fields::AbstractArray{Any}, field_name::Symbol)
    for (i, field) in enumerate(fields)
        if isa(field, Expr) && field.args[1] == field_name 
            return i, field 
        end
    end
    return (nothing, nothing)
end

function get_struct_name(def::Expr)
    if def.args[2] isa Symbol  # Struct that does not subtype abstract type
        def.args[2]
    elseif def.args[2].args[1] isa Symbol # Subtyped struct with no parametric types 
        def.args[2].args[1]
    else
        def.args[2].args[1].args[1]
    end
end

function inject_hash_behaviour(expr, typecode)
    # Assumes last line for the inner constructor is a new

    functions = findall(expr -> expr isa Expr && expr.head === :function, expr.args[3].args)
    for func_i in functions 
        inner_ctor = expr.args[3].args[func_i]

        new_call = inner_ctor.args[2].args[end].args

        # Make sure inner ctor ends l
        (new_call[1] === :new || (new_call[1] isa Expr && new_call[1].args[1] === :new)) || continue

        sig_param_is = findall(param -> param isa Expr && param.head === :macrocall && (param.args[1] === Symbol("@sig") || param.args[1] === Symbol("@sig_t")), new_call)
        sig_exprs = []
        for sig_param_i ∈ sig_param_is
            # sig_param_i will be the index of parameter wrapped in a @sig 
            # we unwrap it and mark that expression needs to be included in the hash
            param = new_call[sig_param_i]

            # sig(Tuple) => hash((t1.hash_code, t2.hash_code, ..., tn.hash_code))
            # sig(expr) => hash(expr) 

            is_tuple = param.args[3] isa Expr && param.args[3].head === :tuple 
            new_call[sig_param_i] = new_call[sig_param_i].args[3]


            if is_tuple 
                with_hashcode_accessors = map(expr -> :(hash_behaviour($expr)), param.args[3].args)
                push!(sig_exprs, Expr(:tuple, with_hashcode_accessors...))
            else 
                push!(sig_exprs, new_call[sig_param_i])
            end 
        end


        hash_expr = :( $(esc(:hash))(($typecode, $(sig_exprs...))) )

        insert!(new_call, 2, hash_expr)
    end
end

"""
    Inserts predefined fields for a struct
    Parses fields tags into metadata
    Injects hash logic
"""
typecode = 0
macro core_type(expr)
    struct_fields = expr.args[3].args

    methods_to_add, fields_to_add = [], []
    struct_name = get_struct_name(expr)

    prepend!(fields_to_add, (field(:hash_code, :UInt32),))

    metadata = Dict{Symbol, Int}()
    tags_index, tags = find_field(struct_fields, :tags)
    if tags_index !== nothing 
        deleteat!(struct_fields, tags_index)

        wanted_tags = tags.args[2].args
        for wanted_tag ∈ wanted_tags 
            metadata[wanted_tag.args[1]] = wanted_tag.args[2]
        end
    end

    # Expand fields variable into respective struct members
    fields_index, fields = find_field(struct_fields, :fields)
    if fields_index !== nothing
        deleteat!(struct_fields, fields_index)

        wanted_fields = fields.args[2].args

        for wanted_field ∈ wanted_fields
            field_name = Symbol(:wanted_field)
    
            field_name ∉ keys(fields) && error("don't recognize field $(field_name)")

            if field_name === :operands 
                num_operands = get!(metadata, :num_ops, -1)

                var_type = num_operands === -1 ? fields[field_name].type : :(NTuple{$(num_operands), AbstractExpr})
                push!(fields_to_add, field(field_name, var_type))
            else
                push!(fields_to_add, field(field_name, fields[field_name].type))
            end
        end
    end 

    tc = global typecode += 1
    push!(methods_to_add, esc(quote typecode(x::$struct_name) = $tc end))
    
    prepend!(struct_fields, fields_to_add)
    
    inject_hash_behaviour(expr, tc)

    return Expr(:block, expr, methods_to_add...)
end

macro attach_hash_operators(e)
    e.head === :struct || error("can only define hash operators on structs")
    struct_name = e.args[2]

    esc(quote
        $e

        Base.hash(x::$struct_name) = (hash ∘ hash_data)(x)
        Base.:(==)(x::$struct_name, y::$struct_name) = hash_data(x) === hash_data(y) 
    end)
end

hash_data(n::Nothing) = nothing

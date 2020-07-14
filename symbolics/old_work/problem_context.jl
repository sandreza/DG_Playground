import Base: +, *, /, -
abstract type change_to_add end
struct test_struct <: change_to_add end

struct test_field
    a
end

+(a::test_field, b::test_field) = a.a + b.a
*(a::test_field, b::test_field) = a.a * b.a
*(a::test_field, b::test_field, ::change_to_add) = a.a + b.a

a = test_field(1)
b = test_field(2)
c = test_field(3)
cta = test_struct()
a+a
a*a
*(a, a, cta)
d = [a b; a b]
dn = [a b b; b a a]
d + d
d * d
d * dn
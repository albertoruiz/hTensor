import Numeric.LinearAlgebra.Multivector

o  = e 4
x = vector[1,0,0,1]
y = o + e 2
z = vector[0,0,1] + e 4

p = x /\ y /\ z

p1 = o /\ x /\ y
p2 = z /\ vector[1,0,1,1] /\ vector[0,1,1,1]

l = o /\ vector[1,1,1,1]

rot = rotor 3 (pi/4) (e 3)

l' = rot * l * rever rot

inh v = v / (v -| e 4) - e 4

main = do
     print l
     print $ l \/ p1
     print $ l \/ p2
     print $ l \/ p
     print $ inh $ l \/ p1
     print $ inh $ l \/ p2
     print $ inh $ l' \/ p1
     print $ inh $ l' \/ p2

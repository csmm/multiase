data = r"""

type: ?
  ! anything	
  template: (>*)

type: ar
  ! Argon atom
  template: (>Ar)

type:c1o
  ! Carbon in CO
  template: [>C[~O]]
end_type

type:c2=
  ! Carbon in =C= (e.g. CO2, CS2)
  template: [>C[~*][~*]]

type:c3'
    ! Carbonyl carbon [one polar substituent such as O,N]
    ! e.g. amide, acid and ester
    template: [>C (~O) (~CH) (~ON)]

type:c3a
  ! SP2 aromatic carbon
  template:[>CR (~*)(~*)(~*)]

type:c4
  ! generic SP3 carbon
  template: [>C(-*)(-*)(-*)(-*)]

type:c41o
  ! Carbon, sp3, in methanol (and dimethyl ether?)
  template: [>C(-O(-CH))(-H)(-H)(-H)]

type: c43
  ! sp3 carbon with 1 h and 3 heavy atoms
  template: (>C(-H)(-^H)(-^H)(-^H))

type:c43o
  ! Carbon, sp3, in secondary alcohols
  template: [>C(-O(-H))(-H)(-C)(-C)]

type: c44
  ! sp3 carbon with four heavy atoms attached
  template: [>C(-^H)(-^H)(-^H)(-^H)]

type: c4o
  ! alpha carbon (e.g. alpha to oxygen in ethers and alcohols)
  template: [>C(-O)(-*)(-*)(-*)]

type: c4z
  ! Carbon, sp3, bonded to -N3 (azides)
  template: [>C(-N(~N(~N)))(-*)(-*)(-*)]

type:h1
  ! nonpolar hydrogen 
  template: (>H (-CSi) )

type:h1h
  ! Hydrogen in H2
  template: [>H[-H]]

type:h1o
  ! strongly polar hydrogen (bonded to fluorine, nitrogen, Oxygen - h* in pcff)
  template: (>H(-ONF))

type: he
  ! Helium atom
  template: (>He)

type: kr
  ! Krypton atom
  template: (>Kr)

type:n1n
  ! Nitrogen in N2
  template: [>N[~N]]

type:n1o
  ! Nitrogen in NO
  template: [>N[~O]]

type:n1z
  ! Nitrogen, terminal atom in -N3
  template: [>N[~N[~N(~*)]]]

type:n2=
  ! Nitrogen (in phosphazenes, or generic???)
  template: [>N(~*)(~*)]

type:n2o
  ! Nitrogen in NO2
  template: [>N[~O][~O]]

type:n2t
  ! Nitrogen, central atom in -N3
  template: [>N[~N][~N(~*)]]

type:n2z
  ! Nitrogen, first atom in -N3
  template: [>N[~N[~N]](~*)]

type: n3m
  ! sp3 nitrogen in amides without hydrogen
  template: [>N(-C[=O])(-C)(-C)]

type: n3o
  ! Nitrogen in nitro group
  template: (>N[~O][~O](-C))

type: ne
  ! Neon atom
  template: (>Ne)

type:o1=
  ! Oxygen in NO2 and SO2 [and carbonyl]
  template: [>O(~NSC)]

type:o1=*
  ! Oxygen in CO2
  template: [>O[~C[~O]]]

type:o12
  ! Oxygen in nitro group -NO2
  template: [>O[~N[~O](~*)]]

type:o1c
  ! Oxygen in CO
  template: [>O[~C]]

type:o1n
  ! Oxygen in NO
  template: [>O[~N]]

type:o1o
  ! Oxygen in O2
  template: [>O[~O]]

type:o2
  ! Generic oxygen with two bonds attached
  template: [>O(~*)(~*)]

# Non-aromatic?
type:o2e
  ! Ether oxygen
  template: [>O(-C)(-C)]

type:o2h
  ! Hydroxyl oxygen
  template: (>O[-H](~*))

type:o2n
  ! Oxygen in nitrates
  template: (>O[~N[~O][~O]](~C))

type:o2s
  ! Ester oxygen
  template: (>O[~C[~O](~*)](~C))
end_type

type: o2z
   ! Oxygen in siloxanes and zeolites
   template: (>O(-Si)(-SiH) )

type: p4=
   ! Phosphorous [in phosphazenes]
   template: (>P(~*)(~*)(~*)(~*))

type:s1=
  ! Sulfur in CS2
  template: [>S[~C[~S]]]

type:s2=
  ! Sulfur in SO2
  template: [>S[~O][~O]]

type: si4
   ! Generic silicon with four bonds attached
   template: (>Si(-*)(-*)(-*)(-*))

type: si4c
   ! A subset of si4, non-hydrogen atom attached [siloxanes??]
   template: (>Si(-O)(-OC)(-OC)(-OC))

type: xe
  ! Xenon atom
  template: (>Xe)
  
 precedence:
(?
  (ar)
  (c1o)
  (c2=)
  (c3a) (c3')
  (c4 (c43 (c43o)) (c44) (c4o(c41o)) (c4z) )
  (h1) (h1h) (h1o)
  (he)
  (kr)
  (n1n) (n1o) (n1z)
  (n2= (n2o) (n2t) (n2z) )
  (n3m) (n3o)
  (ne)
  (o1= (o1=*) (o12) (o1c) (o1n) ) (o1o)
  (o2 (o2e(o2s)) (o2h) (o2n) (o2z) )
  (p4=)
  (s1=)
  (s2=)
  (si4 (si4c) )
  (xe)
)

"""
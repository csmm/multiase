data = """

#type: H0
#	! H aliph. bond. to C with 1 electrwd. group (03GLY)
#	template:

type: Br
	! bromine
	template: (Br)

type: C
	!  sp2 C carbonyl group 
	template: [C [O](*)(*)]
	
type: CA
	! sp2 C pure aromatic (benzene)
	template: [CR6 (*)(C)(C)]
	
type: CB
	! sp2 aromatic C, 5&6 membered ring junction
	template: [CR6 (CR6)(CR5)(CR5)]

#type: CC
#	! sp2 aromatic C, 5 memb. ring HIS
#	template: [CR5 (C)(C)(*)]

#type: CK
#	! sp2 C 5 memb.ring in purines
#	template:

#type: CM
#	! sp2 C  pyrimidines in pos. 5 & 6
#	template:

#type: CN
#	! sp2 C aromatic 5&6 memb.ring junct.(TRP)
#	template:

#type: CQ
#	! sp2 C in 5 mem.ring of purines between 2 N
#	template:

#type: CR
#	! sp2 arom as CQ but in HIS
#	template:

type: CT
	! sp3 aliphatic C
	template: [C (*)(*)(*)(*)]

#type: CV
#	! sp2 arom. 5 memb.ring w/1 N and 1 H (HIS)
#	template:

#type: CW
#	! sp2 arom. 5 memb.ring w/1 N-H and 1 H (HIS)
#	template:

#type: C*
#	! sp2 arom. 5 memb.ring w/1 subst. (TRP)
#	template:

type: C0
	! calcium
	template: (Ca)

type: F
	! fluorine
	template: (F)

type: H
	! H bonded to nitrogen atoms
	template: [H (N)]

type: HC
	! H aliph. bond. to C without electrwd.group
	template: [H [C(CH)(CH)(CH)]]

type: H1
	! H aliph. bond. to C with 1 electrwd. group
	template: [H [C(CH)(CH)(ONSFCl)]]

type: H2
	! H aliph. bond. to C with 2 electrwd.groups
	template: [H [C(CH)(ONSFCl)(ONSFCl)]]

type: H3
	! H aliph. bond. to C with 3 eletrwd.groups
	template: [H [C(ONSFCl)(ONSFCl)(ONSFCl)]]

type: HA
	! H arom. bond. to C without elctrwd. groups
	template: [H [CR(CH)(CH)]]

type: H4
	! H arom. bond. to C with 1 electrwd. group
	template: [H [CR(CH)(ONSFCl)]]

type: H5
	! H arom. bond. to C with 2 electrwd. groups
	template: [H [CR(ONSFCl)(ONSFCl)]]

type: HO
	! hydroxyl group
	template: [H [O(C)]]

type: HS
	! hydrogen bonded to sulphur (pol?)
	template: [H (S)]

type: HW
	! H in TIP3P water
	template: [H [O(H)]]

#type: HP
#	! H bonded to C next to positively charged gr
#	template:

type: I
	! iodine  (Applequist)
	template: (I)

type: Cl
	! chlorine  (Applequist)
	template: (Cl)

type: Na
	! Na+, ions pol:J.PhysC,11,1541,(1978)
	template: (Na)

#type: IB
#	! 'big ion w/ waters' for vacuum (Na+, 6H2O)
#	template:

type: MG
	! magnesium
	template: (Mg)

type: N
	! sp2 nitrogen in amide groups
	template: [N (C[O])(*)(*)]

#type: NA
#	! sp2 N in 5 memb.ring w/H atom (HIS)
#	template:

#type: NB
#	! sp2 N in 5 memb.ring w/LP (HIS,ADE,GUA)
#	template:

#type: NC
#	! sp2 N in 6 memb.ring w/LP (ADE,GUA)
#	template:

type: N2
	! sp2 N in amino groups
	template: [N (*)(*)(*)]

#type: N3
#	! sp3 N for charged amino groups (Lys, etc)
#	template:

#type: N*
#	! sp2 N 
#	template:

type: O
	! carbonyl group oxygen
	template: [O (C)]

type: OW
	! oxygen in TIP3P water
	template: [O (H)(H)]

type: OH
	! oxygen in hydroxyl group
	template: [O (H)(C)]

type: OS
	! ether and ester oxygen
	template: [O (^H)(^H)]

#type: O2
#	! carboxyl and phosphate group oxygen
#	template:

type: P
	! phosphate,pol:JACS,112,8543,90,K.J.Miller
	template: (P)

#type: S
#	! S in disulfide linkage,pol:JPC,102,2399,98
#	template:

#type: SH
#	! S in cystine
#	template:

type: CU
	! copper
	template: (Cu)

type: FE
	! iron
	template: (Fe)

type: K
	! potassium
	template: (K)

type: Rb
	! rubidium
	template: (Rb)

type: Cs
	! cesium
	template: (Cs)

type: Li
	! lithium, ions pol:J.PhysC,11,1541,(1978)
	template: (Li)

type: Zn
	! Zn2+
	template: (Zn)

"""

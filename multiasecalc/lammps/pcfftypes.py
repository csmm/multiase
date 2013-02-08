data = '''
 #2.1 11   Ag     107.86800     Ag          0        Silver metal
 #2.1 11   Al      26.98200     Al          0        Aluminium metal
 #2.1 11   Au     196.96700     Au          0        Gold metal
 #1.0  1   Br      79.90900     Br          1        bromine ion
 #1.0  1   Cl      35.45300     Cl          1        chlorine ion
 #2.1 11   Cr      51.99600     Cr          0        Chromium metal
 #2.1 11   Cu      63.54600     Cu          0        Copper metal
 #2.1 11   Fe      55.84700     Fe          0        Iron metal
 #2.1 11   K       39.10200      K          0        Potassium metal
 #2.1 11   Li       6.94000     Li          0        Lithium metal
 #2.1 11   Mo      95.94000     Mo          0        Molybdenum metal
 #2.1 11   Na      22.99000     Na          0        Sodium metal
 #2.1 11   Ni      58.71000     Ni          0        Nickel metal
 #2.1 11   Pb     207.20000     Pb          0        Lead metal
 #2.1 11   Pd     106.40000     Pd          0        Palladium metal
 #2.1 11   Pt     195.09000     Pt          0        Platinum metal
 #2.1 11   Sn     118.69000     Sn          0        Tin metal
 #2.1 11   W      183.85000      W          0        Tungsten metal
 #2.1  8   ar      39.94400     Ar          0        Argon
 #3.0 10   az      26.98200     Al          4        aluminium atom in zeolites
 #1.0  1   br      79.90900     Br          1        bromine atom
type: c
	! generic SP3 carbon
	template: [C (*)(*)(*)(*)]

#1.0  1   c+      12.01115      C          3        C in guanidinium group
#1.0  1   c-      12.01115      C          3        C in charged carboxylate

type: c1
	! sp3 carbon with 1 H 3 heavies
	template:  [C (H)(^H)(^H)(^H)]
 
type: c2
	! sp3 carbon with 2 H's, 2 Heavy's
	template:  [C (H)(H)(^H)(^H)]
 
type: c3
	! sp3 carbon with 3 hHs 1 heavy
	! or methane
	template:  [C (H)(H)(H)(*)]

type: c3h
	! sp3 carbon in 3-membered ring with hydrogens
	template:  [CR3 (H)(H)(*)(*)]

# 1.0  1   c3m     12.01115      C          4        sp3 carbon in 3-membered ring
# 1.0  1   c4h     12.01115      C          4        sp3 carbon in 4-membered ring with hydrogens
# 1.0  1   c4m     12.01115      C          4        sp3 carbon in 4-membered ring

type: c5
	! sp2 aromatic carbon in 5-membered ring
	template: [CR5 (*)(*)(*)]

type: c=
	! non aromatic end doubly bonded carbon
	template: [C (*)(H)(H)]

type: c=1
	! non aromatic, next to end doubly bonded carbon
	template: [C [C(H)(H)](*)(*)]

type: c=2
	! non aromatic doubly bonded carbon
	template: [C (*)(*)(*)]
 #2.1  8   c_0     12.01115      C          3        carbonyl carbon of aldehydes, ketones
 #2.1  8   c_1     12.01115      C          3        carbonyl carbon of acid, ester, amide 
 #2.1  8   c_2     12.01100      C          3        carbonyl carbon of carbamate, urea
 #1.0  1   c_a     12.01115      C          4        general amino acid alpha carbon (sp3)
 #1.0  1   ca+     40.08000     Ca          1        calcium ion  
 #1.0  1   cg      12.01115      C          4        sp3 alpha carbon in glycine
 #1.0  1   ci      12.01115      C          3        sp2 aromatic carbon in charged imidazole ring (His+)
 #1.0  1   cl      35.45300     Cl          1        chlorine atom
 #1.0  1   co      12.01115      C          4        sp3 carbon in acetals
 #1.0  1   coh     12.01115      C          4        sp3 carbon in acetals with hydrogen
 
type: cp
	! sp2 aromatic carbon
	template: [CR (*)(*)(*)]

 #1.0  1   cr      12.01115      C          3        C in neutral arginine
 #1.0  1   cs      12.01115      C          3        sp2 aromatic carbon in 5 membered ring next to S
 #1.0  1   ct      12.01115      C          2        sp carbon involved in a triple bond
 #2.0  5   cz      12.01100      C          3        carbonyl carbon of carbonate

precedence: (c (c3h) (c1) (c2) (c3))

precedence: ((c5) (cp) (c=2 (c=) (c=1)))
 
 #1.0  1   dw       2.01400      D          1        deuterium in heivy water    
 #2.1  8   f       18.99840      F          1        fluorine  atom  
 
type: h
	! generic hydrogen bound to C, Si,or H    
	template: [H (CHSi)]
	
 #1.0  1   h*       1.00797      H          1        hydrogen bonded to nitrogen, Oxygen
 #1.0  1   h+       1.00797      H          1        charged hydrogen in cations
 #3.0 10   hb       1.00782      H          1        hydrogen atom in bridging hydroxyl group
 #1.0  1   hc       1.00797      H          1        hydrogen bonded to carbon  
 #2.1  8   he       4.00300     He          0        Helium                      
 #1.0  1   hi       1.00797      H          1        Hydrogen in charged imidazole ring
 #1.0  1   hn       1.00797      H          1        hydrogen bonded to nitrogen 
 #2.1  8   hn2      1.00800      H          1        amino hydrogen

type: ho
	! hydrogen bonded to oxygen
	template: [H (O)]


#type: ho2
#	! hydroxyl hydrogen
#	! next to o_2 or oz
#	template

 #3.0 10   hoa      1.00782      H          1        hydrogen atom in terminal hydroxyl group on aluminium
 #3.0 10   hos      1.00782      H          1        hydrogen atom in terminal hydroxyl group on silicon
 #1.0  1   hp       1.00797      H          1        hydrogen bonded to phosphorus 
 #1.0  1   hs       1.00797      H          1        hydrogen bonded to sulfur  
 #2.2  9   hsi      1.00800      H          1        silane hydrogen
 #1.0  1   hw       1.00797      H          1        hydrogen in water
 #1.0  1   i      126.90440      I          1        iodine atom
 #2.1  8   kr      83.80000     Kr          0        Krypton
 #1.0  1   n       14.00670      N          3        generic sp2 nitrogen (in amids))
 #1.0  1   n+      14.00670      N          4        sp3 nitrogen in protonated amines
 #1.0  1   n1      14.00670      N          3        sp2 nitrogen in charged arginine
 #1.0  1   n2      14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2) 
 #1.0  1   n3m     14.00670      N          3        sp3 nitrogen in 3- membered ring
 #1.0  1   n3n     14.00670      N          3        sp2 nitrogen in 3- membered ring
 #1.0  1   n4      14.00670      N          4        sp3 nitrogen in protonated amines
 #1.0  1   n4m     14.00670      N          3        sp3 nitrogen in 4- membered ring
 #1.0  1   n4n     14.00670      N          3        sp2 nitrogen in 4- membered ring
 #1.0  1   n=      14.00670      N          2        non aromatic end doubly bonded nitrogen
 #1.0  1   n=1     14.00670      N          2        non aromatic, next to end doubly bonded carbon
 #1.0  1   n=2     14.00670      N          2        non aromatic doubly bonded nitrogen            
 #1.0  1   n_2     14.01000      N          3        nitrogen of urethane
 #1.0  1   na      14.00670      N          3        sp3 nitrogen in amines
 #1.0  1   nb      14.00670      N          3        sp2 nitrogen in aromatic amines
 #2.1  8   ne      20.18300     Ne          0        Neon 
 #1.0  1   nh      14.00670      N          3        sp2 nitrogen in 5 or 6 membered ring
 #1.0  1   nh+     14.00670      N          3        protonated nitrogen in 6 membered ring
 #1.0  1   nho     14.00670      N          3        sp2 nitrogen in 6 membered ring next to a carbonyl 
 #1.0  1   ni      14.00670      N          3        nitrogen in charged imidazole ring
 #1.0  1   nn      14.00670      N          3        sp2 nitrogen in aromatic amines
 #1.0  1   np      14.00670      N          2        sp2 nitrogen in 5- or 6- membered ring
 #1.0  1   npc     14.00670      N          3        sp2 nitrogen in 5- or 6- membered ring and with a heavy atom
 #1.0  1   nr      14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2)
 #1.0  1   nt      14.00670      N          1        sp nitrogen involved in a triple bond
 #1.0  1   nz      14.00670      N          1        sp3 nitrogen bonded to two atoms

type: o
	! generic SP3 oxygen
	template: [O (*)(*)]
	
type: o*
	! oxygen in water
	template: [O (H)(H)]

 #1.0  1   o-      15.99940      O          1        partial double oxygen  
 #1.0  1   o3e     15.99940      O          2        sp3 oxygen  in three membered ring
 #1.0  1   o4e     15.99940      O          2        sp3 oxygen  in  four  membered ring

type: o=
	! oxygen double bonded to O, C, S, N, P
	template: [O (OCSNP)]

type: o_1
	! oxygen in carbonyl group
	template: [O [C(*)(*)]]

 #2.1  8   o_2     15.99940      O          2        ester oxygen
 #3.0 10   oah     15.99491      O          2        oxygen atom in terminal hydroxyl group on aluminium
 #3.0 10   oas     15.99491      O          2        oxygen atom between aluminium and silicon
 #3.0 10   ob      15.99491      O          3        oxygen atom in bridging hydroxyl group
 #1.0  1   oc      15.99940      O          2        sp3 oxygen  in ether or acetals
 #1.0  1   oe      15.99940      O          2        sp3 oxygen  in ester

type: oh
	! oxygen bonded to hydrogen
	template: [O (H)(*)]

precedence: (o (oh (o*)))

precedence: (o= (o_1))

 #2.0  5   oo      15.99940      O          1        oxygen in carbonyl group, carbonate only
 #1.0  1   op      15.99940      O          2        sp2 aromatic in 5 membered ring 
 #3.0 10   osh     15.99491      O          2        oxygen atom in terminal hydroxyl group on silicon
 #1.0  1   osi     16.00000      O          2        siloxane oxygen
 #3.0 10   oss     15.99491      O          2        oxygen atom betweem two silicons
 #2.0  5   oz      15.99940      O          2        ester oxygen in carbonate
 #1.0  1   p       30.97380      P          4        general phosphorous atom
 #1.0  1   p=      30.97380      P          5        phosphazene phosphorous atom
 #1.0  1   s       32.06400      S          2        sp3 sulfur
 #1.0  1   s'      32.06400      S          1        S in thioketone group
 #1.0  1   s-      32.06400      S          1        partial double sulfur 
 #1.0  1   s1      32.06400      S          2        sp3 sulfur involved in (S-S) group of disulfides
 #1.0  1   s3e     32.06400      S          2        sulfur  in three membered ring
 #1.0  1   s4e     32.06400      S          2        sulfur  in four membered ring
 #1.0  1   sc      32.06400      S          2        sp3 sulfur in methionines (C-S-C) group
 #1.0  1   sf      32.06400      S          1        S in sulfonate group
 #1.0  1   sh      32.06400      S          2        sp3 sulfur in sulfhydryl (-SH) group (e.g. cysteine) 
 #1.0  1   si      28.08600     Si          4        silicon atom
 #1.0  1   sio     28.08600     Si          4        siloxane silicon
 #1.0  1   sp      32.06400      S          2        sulfur in an aromatic ring (e.g. thiophene)
 #3.0 10   sz      28.08600     Si          4        silicon atom in zeolites
 #2.1  8   xe     131.30000     Xe          0        Xenon   
'''
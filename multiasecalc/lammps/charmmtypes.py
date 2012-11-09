# Almost all types for neutral molecules containing C, H, O and N are supported.

data = r"""

# **** Hydrogen ****

type: HGA1
	! aliphatic proton, CH
	template: [H (C(^H)(^H)(^H))]
	
type: HGA2
	! aliphatic proton, CH2
	template: [H (C(H)(^H)(^H))]
	
type: HGA3
	! aliphatic proton, CH3
	! Use also for methane...
	template: [H (C(H)(H)(*))]
	
type: HGA4
	! alkene proton; RHC=
	template: [H [C (C)(CN)]]
	
type: HGA5
	! alkene proton; H2C=CR
	template: [H [C (H)(^H)]]

type: HGA5
	! Use also for HCCH...
	template: [H [C(*)]]

#MASS   261 HGA6     1.00800  ! aliphatic H on fluorinated C, monofluoro
#MASS   262 HGA7     1.00800  ! aliphatic H on fluorinated C, difluoro


# HGAAM 0-2: EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
on NEUTRAL METHYLAMINE groups, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens
type: HGAAM0
	! aliphatic H, NEUTRAL trimethylamine
	template: [H [C [N(C)(C)](*)(*)(*)]]

type: HGAAM1
	! aliphatic H, NEUTRAL dimethylamine
	template: [H [C [N(C)(H)](*)(*)(*)]]

type: HGAAM2
	! aliphatic H, NEUTRAL methylamine
	template: [H [C [N(H)(H)](*)(*)(*)]]

type: HGP1
	! polar H
	template: [H (ON)]

#MASS   267 HGP2     1.00800  ! polar H, +ve charge

type: HGP3
	! polar H, thiol
	template: [H [N(*)]]
	
type: HGP4
	! polar H, neutral conjugated -NH2 group (NA bases)
	template: [H [N(H)[C(*)(*)]]]

#MASS   270 HGP5     1.00800  ! polar H on quarternary ammonium salt (choline)

#HGPAM 1-3: EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
#on NEUTRAL METHYLAMINE groups, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens

type: HGPAM1
	! polar H, NEUTRAL dimethylamine
	template: [H [N(C)(C)]]
	
type: HGPAM2
	! polar H, NEUTRAL methylamine
	template: [H [N(C)(H)]]
	
type: HGPAM3
	! polar H, NEUTRAL ammonia
	template: [H [N(H)(H)]]
	
precedence: (HGP1 (HGP3) (HGP4) (HGPAM1) (HGPAM2) (HGPAM3))

type: HGR51
	! nonpolar H, neutral 5-mem planar ring C, LJ based on benzene
	template: [H [CR5(C)(C)]]
	
type: HGR52
	! Aldehyde H, formamide H (RCOH); nonpolar H, neutral 5-mem planar ring C adjacent to heteroatom or + charge
	template: [H [CR5(^CH)(*)]]

type: HGR52
	! Or carbonyl hydrogen!
	template: [H [C[O](*)]]

#MASS   276 HGR53    1.00800  ! nonpolar H, +ve charge HIS he1(+1)

type: HGR61
	! aromatic H
	template: [H [CR6 (C)(C)]]

type: HGR62
	! nonnpolar H, neutral 6-mem planar ring C adjacent to heteroatom
	template: [>H [CR6 (^C)(*)]]
	
type: HGR62
	! or in a ring containg carbonyl
	template: [H [C R6{(C[O])} (C)(C)]]

precedence: ((HGR61 (HGR62)) (HGR51) (HGR52) (HGA4))

#MASS   279 HGR63    1.00800  ! nonpolar H, NAD+ nicotineamide all ring CH hydrogens
#MASS   280 HGR71    1.00800  ! nonpolar H, neutral 7-mem arom ring, AZUL, azulene, kevo

type: HGTIP3
	! polar H, TIPS3P WATER HYDROGEN
	template: [H [O[H]]]
	
precedence: (HGP1 (HGP3) (HGP4) (HGTIP3))

# **** Carbon ****


type: CG1T1
	! alkyn C
	template: [C (*)(*)]

type: CG1N1
	! C for cyano group
	template: [C [N](*)]

precedence: (CG1T1 (CG1N1))

type: CG2D1
	! alkene; RHC= ; imine C
	template: [C (H)(C)(CN)]

type: CG2D2
	! alkene; H2C=
	template: [>C (H)(H)(*)]


#MASS   285 NG311  12.01100  ! double bond carbon adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC1.
#MASS   286 CG2D2O  12.01100  ! double bond carbon adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC2.
#MASS   287 CG2DC1  12.01100  ! conjugated alkenes, R2C=CR2
#MASS   288 CG2DC2  12.01100  ! conjugated alkenes, R2C=CR2
#MASS   289 CG2DC3  12.01100  ! conjugated alkenes, H2C=
#MASS   290 CG2N1   12.01100  ! conjugated C in guanidine/guanidinium
#MASS   291 CG2N2   12.01100  ! conjugated C in amidinium cation

type: CG2O1
	! carbonyl C: amides
	template: [C [O](N)(CH)]

type: CG2O2
	! carbonyl C: esters, [neutral] carboxylic acids
	template: [C [O](O)(CH)]

type: CG2O4
	! carbonyl C: aldehydes
	template: [C [O](H)(CH)]

type: CG2O5
	! carbonyl C: ketones
	template: [C [O](C)(C)]
	
#MASS   297 CG2O6   12.01100  ! carbonyl C: urea, carbonate
#MASS   298 CG2O7   12.01100  ! CO2 carbon

type: CG2R51
	! 5-mem ring, his CG, CD2(0), trp
	template: [CR5 (*)(*)(*)]

#MASS   300 CG2R52  12.01100  ! 5-mem ring, double bound to N, PYRZ, pyrazole

type: CG2R53
	! 5-mem ring, double bound to N and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo
	template: [>CR5 (-^CH)[-N(-*)](-*)]
	
type: CG2R61
	! 6-mem aromatic C
	template: [CR6 (*)(*)(*)]

type: CG2R62
	! 6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)
	template: [CR6{(C[O])} (*)(*)(*)]

type: CG2R63
	! 6-mem aromatic amide carbon (NA) (and other 6-mem aromatic carbonyls?)
	template: [CR6 [O](N)(CN)]

type: CG2R64
	! 6-mem aromatic amidine and guanidine carbon (between 2 or 3 Ns and double-bound to one of them), NA, PYRM
	template: [CR6 [N(*)](N)(*)]

#MASS   306 CG2R66  12.01100  ! 6-mem aromatic carbon bound to F

#type: CG2R67
#	! 6-mem aromatic carbon of biphenyl
#	template: [CR6 (CR6)(CR6)(CR6)]

type: CG2RC0
	! 6/5-mem ring bridging C, guanine C4,C5, trp
	template: [CR6 (CR5)(CR5)(C)]

# Precedence among sp2 carbons
precedence: (
	(CG2RC0) (CG2R61 (CG2R64) (CG2R63) (CG2R62))
	(CG2R51 (CG2R53)) 
	(CG2D1)
	((CG2O1) (CG2O2) (CG2O4) (CG2O5))
	)

#MASS   309 CG2R71  12.01100  ! 7-mem ring arom C, AZUL, azulene, kevo
#MASS   310 CG2RC7  12.01100  ! sp2 ring connection with single bond(!), AZUL, azulene, kevo

#MASS   311 CG301   12.01100  ! aliphatic C, no hydrogens, neopentane
#MASS   312 CG302   12.01100  ! aliphatic C, no hydrogens, trifluoromethyl
	
type: CG311
	! aliphatic C with 1 H, CH
	template: [C (^H)(^H)(^H)[H]]

#MASS   314 CG312   12.01100  ! aliphatic C with 1 H, difluoromethyl
#MASS   315 CG314   12.01100  ! aliphatic C with 1 H, adjacent to positive N (PROT NTER) (+)

type: CG321
	! aliphatic C for CH2
	template: [>C (^H)(^H)[H][H]]

#MASS   317 CG322   12.01100  ! aliphatic C for CH2, monofluoromethyl
#MASS   318 CG323   12.01100  ! aliphatic C for CH2, thiolate carbon
#MASS   319 CG324   12.01100  ! aliphatic C for CH2, adjacent to positive N (piperidine) (+)

type: CG331
	! aliphatic C for methyl group (-CH3)
	template: [C (*)[H][H][H]] 

#MASS   321 CG334   12.01100  ! aliphatic C for methyl group (-CH3), adjacent to positive N (PROT NTER) (+)

#CG3AM 0-2: EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
#on NEUTRAL METHYLAMINE groups, NOT ETHYL, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens

type: CG3AM0
	! aliphatic C for CH3, NEUTRAL trimethylamine methyl carbon
	template: [C (H)(H)(H)[N(C)(C)]]
	
type: CG3AM1
	! aliphatic C for CH3, NEUTRAL dimethylamine methyl carbon
	template: [C (H)(H)(H)[N(C)(^C)]]
	
type: CG3AM2
	! aliphatic C for CH3, NEUTRAL methylamine methyl carbon
	template: [C (H)(H)(H)[N(^C)(^C)]]

#MASS   325 CG3C31  12.01100  ! cyclopropyl carbon
#!MASS   326 CG3C41  12.01100  ! cyclobutyl carbon RESERVED!
#MASS   327 CG3C50  12.01100  ! 5-mem ring aliphatic quaternary C (cholesterol, bile acids)
#MASS   328 CG3C51  12.01100  ! 5-mem ring aliphatic CH  (proline CA, furanoses)
#MASS   329 CG3C52  12.01100  ! 5-mem ring aliphatic CH2 (proline CB/CG/CD, THF, deoxyribose)
#MASS   330 CG3C53  12.01100  ! 5-mem ring aliphatic CH  adjacent to positive N (proline.H+ CA) (+)
#MASS   331 CG3C54  12.01100  ! 5-mem ring aliphatic CH2 adjacent to positive N (proline.H+ CD) (+)
#MASS   332 CG3RC1  12.01100  ! bridgehead in bicyclic systems containing at least one 5-membered or smaller ring
#!(+) Includes protonated Shiff base (NG3D5, NG2R52 in 2HPP) but NOT amidinium (NG2R52 in IMIM), guanidinium


# **** Nitrogen ****

type: NG1T1
	! N for cyano group
	template: [N(C)]
	
type: NG2D1
	! N for neutral imine/Schiff's base (C=N-R, acyclic amidine, gunaidine)
	template: [N(C)(*)]

#type: NG2S0
#	! N,N-disubstituted amide, proline N (CO=NRR')
	
type: NG2S1
	! peptide nitrogen (CO=NHR)
	template: [N (H)(C)(C[O])]
	
type: NG2S2
	! terminal amide nitrogen (CO=NH2)
	template: [N (H)(H)(C[O])]

type: NG2S3
	! external amine ring nitrogen (planar/aniline), phosphoramidate
	template: [N (H)(H)[CR (*)(*)]]
	
type: NG2O1
	! NITB, nitrobenzene
	template: [N [O][O](CR6)]
	
#MASS   340 NG2P1   14.00700  ! N for protonated imine/Schiff's base (C=N(+)H-R, acyclic amidinium, guanidinium)

type: NG2R50
	! double bound neutral 5-mem planar ring, purine N7
	template: [NR5 (*)(*)]

# This is a tough one: all atoms in ring sp2...
type: NG2R51
	! single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)
	template: [N R5{[C(*)(*)(*)][C(*)(*)(*)][C(*)(*)(*)][C(*)(*)(*)]} (C)(C)(CH)]

type: NG2R51
	! single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)
	template: [N R5{[C(*)(*)(*)][C(*)(*)(*)][C(*)(*)(*)][N(C)(C)]}  (C)(C)(CH)]

precedence: ((NG2R51) (NG3C51))

#MASS   343 NG2R52  14.00700  ! protonated schiff base, amidinium, guanidinium in 5-membered ring, HIS, 2HPP, kevo
#MASS   344 NG2R53  14.00700  ! amide in 5-memebered NON-SP2 ring (slightly pyramidized), 2PDO, kevo

type: NG2R60
	! double bound neutral 6-mem planar ring, pyr1, pyzn
	template: [NR6 (*)(*)]
	
type: NG2R61
	! single bound neutral 6-mem planar ring imino nitrogen; glycosyl linkage
	! Actually seems to only used for ring N next to carbonyl...
	template: [NR6 (C[O])(C)(*)]

precedence: ((NG2R61) (NG2S1))

type: NG2R62
	! double bound 6-mem planar ring with heteroatoms in o or m, pyrd, pyrm
	template: [NR6{(^C)} (*)(*)]
	
type: NG2RC0
	! 6/5-mem ring bridging N, indolizine, INDZ, kevo
	template: [NR6 (CR5)(CR5)(C)]

# sp2 Nitrogen precedence
precedence: (
	(NG2R60 (NG2R62) (NG2RC0))
	(NG2R50)
	(NG2D1)
	)

type: NG301
	! neutral trimethylamine nitrogen
	template: [N (C(*)(*)(*))(C(*)(*)(*))(C(*)(*)(*))]
	
type: NG311
	! neutral dimethylamine nitrogen
	template: [N (C(*)(*)(*))(C(*)(*)(*))(H)]
	
type: NG321
	! neutral methylamine nitrogen
	template: [N (C(*)(*)(*))(H)(H)]
	
type: NG331
	! neutral ammonia nitrogen
	template: [N (H)(H)(H)]
	
type: NG3C51
	! secondary sp3 amine in 5-membered ring
	template: (NR5 (H))

type: NG3N1
	! N in hydrazine, HDZN
	template: [N (N(H))(H)(*)]
	
#MASS   355 NG3P0   14.00700  ! quarternary N+, choline
#MASS   356 NG3P1   14.00700  ! tertiary NH+ (PIP)
#MASS   357 NG3P2   14.00700  ! secondary NH2+ (proline)
#MASS   358 NG3P3   14.00700  ! primary NH3+, phosphatidylethanolamine

# **** Oxygen ****

type: OG2D1
	! carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, uera
	template: [O(C(^C))]
	
#MASS   360 OG2D2   15.99940  ! carbonyl O: negative groups: carboxylates, carbonate

type: OG2D3
	! carbonyl O: ketones
	template: [O[C(C)(C)]]

type: OG2D4
	! 6-mem aromatic carbonyl oxygen (nucleic bases)
	template: [O(CR6)]
	
precedence: ((OG2D4) (OG2D1) (OG2D3))
	
type: OG2D5
	! CO2 oxygen
	template: [O [C[O]]]

#MASS   364 OG2N1   15.99940  ! NITB, nitrobenzene
#MASS   365 OG2P1   15.99940  ! =O in phosphate or sulfate

type: OG2R50
	! FURA, furan
	template: (OR5)
	
#MASS   367 OG3R60  15.99940  ! O in 6-mem cyclic enol ether (PY01, PY02) or ester

type: OG301
	! ether -O-
	template: [O (C)(C)]

type: OG302
	! ester -O-
	template: [O (C[O])(C)]
	
precedence: (OG301 (OG302))

#MASS   370 OG303   15.99940  ! phosphate/sulfate ester oxygen
#MASS   371 OG304   15.99940  ! linkage oxygen in pyrophosphate/pyrosulphate

type: OG311
	! hydroxyl oxygen
	template: [O(H)(^H)]

#MASS   373 OG312   15.99940  ! ionized alcohol oxygen
#MASS   374 OG3C51  15.99940  ! 5-mem furanose ring oxygen (ether)
#MASS   375 OG3C61  15.99940  ! DIOX, dioxane, ether in 6-membered ring !SHOULD WE MERGE THIS WITH OG3R60???

type: OGTIP3
	! TIPS3P WATER OXYGEN
	template: [O(H)(H)]

"""
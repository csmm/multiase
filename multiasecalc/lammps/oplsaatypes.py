data = """

type: opls_111
	! O TIP3P Water
	template: [O (H)(H)]
	
type: opls_112
	! H TIP3P Water
	template: [H [O(H)]]

type: opls_135
	! alkane CH3
	template: [C (C)(H)(H)(H)]
	
type: opls_136
	! alkane CH2 
	template: [C (C)(C)(H)(H)]
	
type: opls_137
	! alkane CH 
	template: [C (C)(C)(C)(H)]
	
type: opls_138
	! alkane CH4 
	template: [C (H)(H)(H)(H)]
	
# opls_139   12.01100  ; alkane C   

type: opls_140
	! alkane H
	template: [H [C (*)(CH)(CH)]]
	
type: opls_141
	! alkene C (R2-C=) 
	template: [C (^H)(^H)(C)]
	
type: opls_142
	! alkene C (RH-C=) 
	template: [C (^H)(^H)(H)]
	
 type: opls_143
	! alkene C (H2-C=) 
	template: [C (^H)(H)(H)]
	
type: opls_144
	! alkene H (H-C=)
	template: [H [C(*)(*)]]

type: opls_145
	! Benzene C - 12 site JACS,112,4768-90
	template: [CR6 (C)(C)(H)]

# This type doesn't make sense, it should not have a charge...
# opls_145B  12.01100  ; Biphenyl C1

type: opls_146
	! Benzene H - 12 site.
	template: [H [CR6 (C)(C)]]
	
type: opls_147
	! Naphthalene fusion C (C9)
	template: [CR6 (CR6)(CR6)(CR6)]

type: opls_148
	! C: CH3, toluene
	template: [C (CR6)(H)(H)(H)]

# hydrogens
precedence: ((opls_146) (opls_144))

type: opls_149
	! C: CH2, ethyl benzene
	template: [C (CR6)[C(H)(H)(H)](H)]

type: opls_150
	! diene =CH-CH=; use #178 for =CR-CR=
	template: [C (H)[C(*)(*)][C(C)(H)]]
	
# opls_151   35.45300  ; Cl in alkyl chlorides
# opls_152   12.01100  ; RCH2Cl in alkyl chlorides
# opls_153    1.00800  ; H in RCH2Cl in alkyl chlorides

type: opls_154 
	! all-atom O: mono alcohols
	template: [O (H)(C)]
	
type: opls_155
	! all-atom H(O): mono alcohols, OP(=O)2
	template: [H [O(C)]]
	
type: opls_156
	! all-atom H(C): methanol
	template: [H [C[O(H)](H)(H)]]

type: opls_157
	! all-atom C: CH3 & CH2, alcohols
	template: [C [O(H)](CH)(H)(H)]
	
type: opls_158
	! all-atom C: CH, alcohols
	template: [C [O(H)](C)(C)(H)]
	
type: opls_159
	! all-atom C: C, alcohols
	template: [C [O(H)](C)(C)(C)]
	
#opls_160   12.01100  ; CH2 Trifluoroethanol
#opls_161   12.01100  ; CF3 Trifluoroethanol
#opls_162   15.99940  ; OH  Trifluoroethanol
#opls_163    1.00800  ; HO  Trifluoroethanol
#opls_164   18.99840  ; F   Trifluoroethanol
#opls_165    1.00800  ; H   Trifluoroethanol

type: opls_166
	! C(OH) phenol  Use with all
	template: [CR6 [O(H)](C)(C)]

	
# sp2 carbons
precedence: (
	(opls_166) (opls_147) (opls_145) # Aromatics
	(opls_178) (opls_150)            # Dienes
	(opls_141) (opls_142)            # Alkanes
	)
	
type: opls_167
	!O     phenol
	template: [O (H)(CR6)]

# hydroxyl oxygens
precedence: ((opls_167) (opls_154))

# hydroxyl hydrogens
precedence: ((opls_168) (opls_155))
	
type: opls_168
	! H     phenol  
	template: [H [O(CR6)]]

# *** It is impossible to detect arbitrary polyols with the template syntax ***
# opls_169   15.99940  ; O:    diols 
# opls_170    1.00800  ; H(O): diols
# opls_171   15.99940  ; O:    triols
# opls_172    1.00800  ; H(O): triols
# opls_173   12.01100  ; C(H2OH): triols
# opls_174   12.01100  ; C(HROH): triols
# opls_175   12.01100  ; C(R2OH): triols
# opls_176    1.00800  ; H(CXOH): triols
 
type: opls_178
	! diene =CR-CR=; use #150 for =CH-CH=
	template: [C (C)[C(*)(*)][C(C)(C)]]

# opls_179   15.99940  ; O: anisole
 
type: opls_180
	! O: dialkyl ether
	template: [O (C)(C)]

type: opls_181
	! C(H3OR): methyl ether
	template: [C (O(C))(H)(H)(H)]
	
type: opls_182
	! C(H2OR): ethyl  ether
	template: [C (O(C))(C)(H)(H)]

type: opls_183
	! C(HOR):  i-Pr   ether, allose
	template: [C (O(C))(C)(C)(H)]
	
type: opls_184
	! C(OR):   t-Bu   ether
	template: (C (O(C))(C)(C)(C))
 opls_185    1.00800  ; H(COR): alpha H ether
 
 """
 
"""
 opls_186   15.99940  ; O: acetal ether 
 opls_187   15.99940  ; O(H): hemiacetal
 opls_188    1.00800  ; H(O): hemiacetal
 opls_189   12.01100  ; C(H2O2): acetal OCH2O
 opls_190    1.00800  ; H(CHO2): acetal OCH2O
 opls_191   12.01100  ; C(H2O2): hemiacetal OCH2OH
 opls_192    1.00800  ; H(CHO2): hemiacetal OCH2OH
 opls_193   12.01100  ; C(HCO2): acetal OCHRO
 opls_194    1.00800  ; H(CHO2): acetal OCHRO
 opls_195   12.01100  ; C(HCO2): hemiacetal OCHROH
 opls_196    1.00800  ; H(C2O2): hemiacetal OCHROH
 opls_197   12.01100  ; C(C2O2): acetal OCRRO
 opls_198   12.01100  ; C(C2O2): hemiacetal OCRROH
 opls_199   12.01100  ; C(O,Me): anisole
 opls_200   32.06000  ; all-atom S: thiols
 opls_201   32.06000  ; S   IN H2S  JPC,90,6379 (1986)
 opls_202   32.06000  ; all-atom S: sulfides, S=C
 opls_203   32.06000  ; all-atom S: disulfides
 opls_204    1.00800  ; all-atom H(S): thiols
 opls_205    1.00800  ; H   IN H2S  JPC,90,6379 (1986)
 opls_206   12.01100  ; all-atom C: CH2, thiols
 opls_207   12.01100  ; all-atom C: CH, thiols
 opls_208   12.01100  ; all-atom C: C, thiols
 opls_209   12.01100  ; all-atom C: CH3, sulfides
 opls_210   12.01100  ; all-atom C: CH2, sulfides
 opls_211   12.01100  ; all-atom C: CH, sulfides
 opls_212   12.01100  ; all-atom C: C, sulfides
 opls_213   12.01100  ; all-atom C: CH3, disulfides
 opls_214   12.01100  ; all-atom C: CH2, disulfides
 opls_215   12.01100  ; all-atom C: CH, disulfides
 opls_216   12.01100  ; all-atom C: C, disulfides
 opls_217   12.01100  ; all-atom C: CH3, methanethiol
 opls_218   12.01100  ; C in CH2OH - benzyl alcohols
 opls_219   12.01100  ; C in CHROH - benzyl alcohols
 opls_220   12.01100  ; C in CR2OH - benzyl alcohols
 opls_221   12.01100  ; C(CH2OH)   - benzyl alcohols
 opls_222   32.06000  ; S in thioanisoles
 opls_223   12.01100  ; C in RCH2NH2. Use #223B for AA Calpha.
 opls_223B  12.01100  ; Gly Calpha
 opls_224   12.01100  ; C in R2CHNH2. Use #224B for AA Calpha.
 opls_224B  12.01100  ; Calpha in most AA (except Gly,Pro,Aib)
 opls_225   12.01100  ; C in R3CNH2. Use #225B for AA Calpha.
 opls_225B  12.01100  ; Aib Calpha.
 opls_226   35.45300  ; chloroalkene Cl (ClH-C=) - see also #398
 opls_227   12.01100  ; chloroalkene C (ClH-C=)
 opls_228   12.01100  ; C(SMe)  thioanisole
 opls_229   12.01100  ; C on N: secondary N-CHR2 amide
 opls_230   12.01100  ; C on N: secondary N-CR3  amide
 opls_231   12.01100  ; C: C=O in benzophenone
 opls_232   12.01100  ; C: C=O in benzaldehyde,acetophenone (CH)
 opls_233   12.01100  ; C: C=O in acetophenone (CMe)
 opls_234   12.01100  ; C: C=O in benzamide
 opls_235   12.01100  ; C=O in amide, dmf, peptide bond
 opls_236   15.99940  ; O: C=O in amide.   Acyl R on C in amide is neutral - 
 opls_237   14.00670  ; N: primary amide.  use alkane parameters.
 opls_238   14.00670  ; N: secondary amide, peptide bond (see #279 for formyl H)
 opls_239   14.00670  ; N: tertiary amide
 opls_240    1.00800  ; H on N: primary amide
 opls_241    1.00800  ; H on N: secondary amide
 opls_242   12.01100  ; C on N: secondary N-Me amide
 opls_243   12.01100  ; C on N: tertiary  N-Me amide
 opls_244   12.01100  ; C on N: secondary N-CH2R amide
 opls_245   12.01100  ; C on N: tertiary  N-CH2R amide, Pro CD
 opls_246   12.01100  ; C on N: tertiary  N-CHR2 amide, Pro CA
 opls_247   12.01100  ; C in O=C(NH2)2  Urea
 opls_248   15.99940  ; O in O=C(NH2)2  Urea Isr. J. Chem
 opls_249   14.00670  ; N in O=C(NH2)2  Urea 33, 323 (93)
 opls_250    1.00800  ; H in O=C(NH2)2  Urea
 opls_251   14.00670  ; N   in imide
 opls_252   12.01100  ; C(=O) in imide
 opls_253   15.99940  ; O   in imide
 opls_254    1.00800  ; H(N) in imide
 opls_255    1.00800  ; H(C) in formimide
 opls_256   12.01100  ; C in CH3  imide        
 opls_257   12.01100  ; C in RCH2 imide     
 opls_258   12.01100  ; C in R2CH imide
 opls_259   12.01100  ; C in R3C  imide
 opls_260   12.01100  ; C(CN)  benzonitrile
 opls_261   12.01100  ; C(N)   benzonitrile
 opls_262   14.00670  ; N      benzonitrile
 opls_263   12.01100  ; C(Cl)  chlorobenzene
 opls_264   35.45300  ; Cl     chlorobenzene
 opls_265   14.00670  ; N: N-phenylacetamide
 opls_266   12.01100  ; ipso C in N-phenylacetamide
 opls_267   12.01100  ; Co in CCOOH carboxylic acid 
 opls_268   15.99940  ; Oh in CCOOH R in RCOOH is
 opls_269   15.99940  ; Oc in CCOOH neutral; use #135-#140
 opls_270    1.00800  ; H  in CCOOH
 opls_271   12.01100  ; C in COO- carboxylate
 opls_272   15.99940  ; O: O in COO- carboxylate,peptide terminus
 opls_273   12.01100  ; C: CH3, carboxylate ion
 opls_274   12.01100  ; C: CH2, carboxylate ion
 opls_275   12.01100  ; C: CH,  carboxylate ion
 opls_276   12.01100  ; C: C,   carboxylate ion
 opls_277   12.01100  ; AA C: aldehyde - for C-alpha use #135-#139
 opls_278   15.99940  ; AA O: aldehyde  
 opls_279    1.00800  ; AA H-alpha in aldehyde & formamide
 opls_280   12.01100  ; AA C: ketone - for C-alpha use #135-#139
 opls_281   15.99940  ; AA O: ketone 
 opls_282    1.00800  ; AA H on C-alpha in ketone & aldehyde
 opls_283   12.01100  ; CA on C-terminal ALA,CYS,SER,THR,HIS,ASP,ASN
 opls_284   12.01100  ; CA on C-terminal GLY
 opls_285   12.01100  ; CA on C-terminal PRO
 opls_286   14.00670  ; N (NH4+) JPC,90,2174 (1986)
 opls_287   14.00670  ; N (RNH3+) JPC,90,2174 (1986)
 opls_288   14.00670  ; N (R4N+)  JPC,90,2174 (1986)  
 opls_289    1.00800  ; H (NH4+)  JPC,90,2174 (1986)  
 opls_290    1.00800  ; H (RNH3+) JPC,90,2174 (1986)
 opls_291   12.01100  ; C in  CH3NH3+      
 opls_292   12.01100  ; C in  RCH2NH3+
 opls_292B  12.01100  ; CA in GLY-NH3+ N-term.
 opls_293   12.01100  ; C in  R2CHNH3+     
 opls_293B  12.01100  ; CA in NH3+ N-term, All AA except GLY & PRO
 opls_294   12.01100  ; C in  R3CNH3+      
 opls_295   12.01100  ; AA C-alpha on N-term PRO
 opls_296   12.01100  ; AA:C-delta in N-term PRO NH2+
 opls_297   12.01100  ; CT in  CH3NH2+R
 opls_298   12.01100  ; AA C-alpha in Gly zwitterion
 opls_299   12.01100  ; AA C-alpha in Ala zwitterion
 opls_300   14.00670  ; N: guanidinium NH2
 opls_301    1.00800  ; H: guanidinium NH2
 opls_302   12.01100  ; C: guanidinium C+
 opls_303   14.00670  ; N: guanidinium NHR
 opls_304    1.00800  ; H: guanidinium NHR
 opls_305   12.01100  ; C: CH3, methylguanidinium  
 opls_306   12.01100  ; C: CH3, ethylguanidinium
 opls_307   12.01100  ; C: CH2(D), ARG, ethylguanidinium
 opls_308   12.01100  ; C: CH2(G), ARG
 opls_309   14.00670  ; N (R2NH2+), N-terminal PRO NH2+
 opls_310    1.00800  ; H (R2NH2+)
 opls_311   14.00670  ; DAP N1   (Diaminopyridine)
 opls_312   12.01100  ; DAP C2   
 opls_313   14.00670  ; DAP N-amine
 opls_314    1.00800  ; DAP H-amine 
 opls_315   12.01100  ; DAP C3
 opls_316    1.00800  ; DAP H3
 opls_317   12.01100  ; DAP C4
 opls_318    1.00800  ; DAP H4
 opls_319   14.00670  ; Uracil & Thymine N1 - use #319B for nucleoside
 opls_319B  14.00670  ; Uracil & Thymine N1 - only for nucleoside
 opls_320   12.01100  ; Uracil & Thymine C2
 opls_321   14.00670  ; Uracil & Thymine N3
 opls_322   12.01100  ; Uracil & Thymine C4
 opls_323   12.01100  ; Uracil & Thymine C5
 opls_324   12.01100  ; Uracil & Thymine C6
 opls_325    1.00800  ; Uracil & Thymine H-N1
 opls_326   15.99940  ; Uracil O-C2
 opls_327    1.00800  ; Uracil H-N3
 opls_328   15.99940  ; Uracil O-C4
 opls_329    1.00800  ; Uracil H-C5
 opls_330    1.00800  ; Uracil H-C6
 opls_331   12.01100  ; Thymine C-C5
 opls_332    1.00800  ; Thymine H-CC5
 opls_333   14.00670  ; Cytosine N1 -use #333B for nucleoside
 opls_333B  14.00670  ; Cytosine N1 - for nucleoside
 opls_334   12.01100  ; Cytosine C2
 opls_335   14.00670  ; Cytosine N3
 opls_336   12.01100  ; Cytosine C4     Nucleotide base
 opls_337   12.01100  ; Cytosine C5     parameters:
 opls_338   12.01100  ; Cytosine C6     JACS,113,2810(1991)
 opls_339    1.00800  ; Cytosine H-N1
 opls_340   15.99940  ; Cytosine O-C2
 opls_341   14.00670  ; Cytosine N-C4
 opls_342    1.00800  ; Cytosine H-NC4/N3
 opls_343    1.00800  ; Cytosine H-NC4/C5
 opls_344    1.00800  ; Cytosine H-C5
 opls_345    1.00800  ; Cytosine H-C6
 opls_346   14.00670  ; Adenine N1
 opls_347   12.01100  ; Adenine C2
 opls_348   14.00670  ; Adenine N3
 opls_349   12.01100  ; Adenine C4
 opls_350   12.01100  ; Adenine C5
 opls_351   12.01100  ; Adenine C6
 opls_352   14.00670  ; Adenine & Guanine N7 
 opls_353   12.01100  ; Adenine & Guanine C8 
 opls_354   14.00670  ; Adenine & Guanine N9 - use #354B for nucleoside
 opls_354B  14.00670  ; Adenine & Guanine N9 - nucleoside only
 opls_355    1.00800  ; Adenine & Guanine H-C2        
 opls_356   14.00670  ; Adenine & Guanine N-C6
 opls_357    1.00800  ; Adenine & Guanine H-NC6/N1
 opls_358    1.00800  ; Adenine & Guanine H-NC6/C5
 opls_359    1.00800  ; Adenine & Guanine H-C8 Guanine
 opls_360    1.00800  ; Adenine & Guanine H-N9 Guanine
 opls_361   14.00670  ; Guanine N1
 opls_362   12.01100  ; Guanine C2
 opls_363   14.00670  ; Guanine N3
 opls_364   12.01100  ; Guanine C4
 opls_365   12.01100  ; Guanine C5
 opls_366   12.01100  ; Guanine C6
 opls_367    1.00800  ; Guanine H-N1
 opls_368   14.00670  ; Guanine N-C2
 opls_369    1.00800  ; Guanine H-NC2
 opls_370   15.99940  ; Guanine O-C6
 opls_371   12.01100  ; 9-Me Adenine or Guanine C-N9
 opls_372    1.00800  ; 9-Me Adenine or Guanine H-CN9
 opls_373   12.01100  ; 1-Me Uracil or Thymine C-N1
 opls_374    1.00800  ; 1-Me Uracil or Thymine H-CN1
 opls_375   12.01100  ; 1-Me Cytosine C-N1
 opls_376    1.00800  ; 1-Me Cytosine H-CN1
 opls_377   14.00670  ; CytH+ N1 Use #377B for nucleoside.
 opls_377B  14.00670  ; CytH+ N1 - nucleoside only
 opls_378   12.01100  ; CytH+ C2      
 opls_379   14.00670  ; CytH+ N3 Protonated cytosine.
 opls_380   12.01100  ; CytH+ C4
 opls_381   12.01100  ; CytH+ C5
 opls_382   12.01100  ; CytH+ C6
 opls_383    1.00800  ; CytH+ H-N1
 opls_384   15.99940  ; CytH+ O-C2
 opls_385    1.00800  ; CytH+ H-N3
 opls_386   14.00670  ; CytH+ N-C4
 opls_387    1.00800  ; CytH+ H-NC4/N3
 opls_388    1.00800  ; CytH+ H-NC4/C5
 opls_389    1.00800  ; CytH+ H-C5
 opls_390    1.00800  ; CytH+ H-C6
 opls_391   12.01100  ; 1-Me CytH+ C-N1
 opls_392    1.00800  ; 1-Me CytH+ H-CN1
 opls_393   30.97376  ; P    dimethylphosphate anion  UA - see #440 for AA
 opls_394   15.99940  ; O(=) dimethylphosphate anion  UA - see #440 for AA
 opls_395   15.99940  ; O(-) dimethylphosphate anion  UA - see #440 for AA
 opls_396   12.01100  ; C in CH3 dimethylphosphate anion UA - see #440 for AA
 opls_400   18.99840  ; F-  JACS 106, 903 (1984)
 opls_401   35.45300  ; Cl- JACS 106, 903 (1984)
 opls_402   79.90400  ; Br- JACS 107, 7793(1985)
 opls_403  126.90450  ; I-  JACS 120, 5104(1998)
 opls_404    6.94100  ; Li+ JACS 106, 903 (1984)
 opls_405   22.98977  ; Na+ JACS 106, 903 (1984)
 opls_406    6.94100  ; Li+
 opls_407   22.98977  ; Na+    Aqvists cation
 opls_408   39.09830  ; K+     parameters:
 opls_409   85.46780  ; Rb+    JPC,94, 8021 (90)
 opls_410  132.90540  ; Cs+
 opls_411   24.30500  ; Mg++
 opls_412   40.08000  ; Ca++
 opls_413   87.62000  ; Sr++
 opls_414  137.33000  ; Ba++
 opls_415   12.01100  ; C  in  CH3S-  thiolate
 opls_416    1.00800  ; H  in  CH3S-
 opls_417   32.06000  ; S  in  CH3S-
 opls_418   12.01100  ; C  in  CH3O-  alkoxide
 opls_419    1.00800  ; H  in  CH3O-
 opls_420   15.99940  ; O  in  CH3O-
 opls_421   12.01100  ; C1 in  CH2CN-  RCN-
 opls_422    1.00800  ; H  in  CH2CN-
 opls_423   12.01100  ; C2 in  CH2CN-   JACS 111,4190 (89)
 opls_424   14.00670  ; N  in  CH2CN-   
 opls_425   12.01100  ; C  in  CH3NH-
 opls_426    1.00800  ; HC in  CH3NH-  RNH-
 opls_427   14.00670  ; N  in  CH3NH-
 opls_428    1.00800  ; HN in  CH3NH-
 opls_429   12.01100  ; C2 in  CH3CH2- RCH2-
 opls_430    1.00800  ; H  in  CH3CH2-
 opls_431   12.01100  ; C1 in  CH3CH2-
 opls_432    1.00800  ; H1 in  CH3CH2-
 opls_433    0.00000  ; LP in  CH3CH2-
 opls_434   15.99940  ; O in OH-  Hyroxide O-H = 0.953 A
 opls_435    1.00800  ; H in OH-  JACS 108, 2517 (86)
 opls_436    0.00000  ; U in UO2+ J Mol Struct 366, 55 (96)
 opls_437   15.99940  ; O in UO2+ r(U-O) = 1.80 A
 opls_440   30.97376  ; P in  Me2PO4-, Me2PO4H
 opls_441   15.99940  ; O= in Me2PO4-, Me2PO4H
 opls_442   15.99940  ; OMe in Me2PO4-, Me2PO4H   dimethylphosphate
 opls_443   12.01100  ; C  in Me2PO4-, Me2PO4H   dimetylphosphate
 opls_444    1.00800  ; H  in Me2PO4-, Me2PO4H    6-31+G* CHELPG
 opls_445   30.97376  ; P  in MeOPO3--, MeOPO3H2
 opls_446   15.99940  ; O= in MeOPO3--, MeOPO3H2
 opls_447   15.99940  ; OMe in MeOPO3--, MeOPO3H2  methyl phosphate
 opls_448   12.01100  ; C  in MeOPO3--, MeOPO3H2   6-31+G* CHELPG
 opls_449    1.00800  ; H  in MeOPO3--, MeOPO3H2
 opls_450   30.97376  ; P  in MePO3Me-, MePO3HMe
 opls_451   15.99940  ; O= in MePO3Me-, MePO3HMe
 opls_452   15.99940  ; OMe in MePO3Me-, MePO3HMe     methyl
 opls_453   12.01100  ; C(O) MePO3Me-, MePO3HMe     methylphosphonate
 opls_454    1.00800  ; H(CO) MePO3Me-, MePO3HMe     6-31+G* CHELPG
 opls_455   12.01100  ; C(P) MePO3Me-, MePO3HMe
 opls_456    1.00800  ; H(CP) MePO3Me-, MePO3HMe
 opls_457   12.01100  ; Cipso  benzyl methylphosphonate
 opls_458   12.01100  ; C(O) benzyl methylphosphonate 
 opls_459    1.00800  ; H(CO) benzyl methylphosphonate
 opls_460   12.01100  ; Cipso  methyl benzylphosphonate
 opls_461   12.01100  ; C(P)  methyl benzylphosphonate
 opls_462    1.00800  ; H(CP) methyl benzylphosphonate
 opls_463   12.01100  ; Cipso C6H5OPO3(2-)  use with #445-#447
 opls_465   12.01100  ; AA C:   esters - for R on C=O, use #280-#282
 opls_466   15.99940  ; AA =O:  esters   
 opls_467   15.99940  ; AA -OR: ester 
 opls_468   12.01100  ; methoxy C in esters - see also #490-#492
 opls_469    1.00800  ; methoxy Hs in esters
 opls_470   12.01100  ; Co in benzoic acid
 opls_471   12.01100  ; Co in methyl benzoate, aryl ester
 opls_472   12.01100  ; Cipso phenyl ester
 opls_473   15.99940  ; AA -OR phenyl ester
 opls_474   32.06000  ; S in sulfonamide, S(=O)2(OR)
 opls_475   15.99940  ; O in sulfonamide, S(=O)2(OR)
 opls_476   12.01100  ; CH3 attached to S of sulfonamide
 opls_477    1.00800  ; H of Me attached to S of sulfonamide
 opls_478   14.00670  ; N: primary amide of sulfonamide
 opls_479    1.00800  ; H on N: primary sulfonamide
 opls_480   14.00670  ; N secondary amide of sulfonamide
 opls_481    1.00800  ; H on N: secondary sulfonamide
 opls_482   12.01100  ; alpha CH3-N of sulfonamide
 opls_483    1.00800  ; H of alpha CH3-N of sulfonamide
 opls_484   12.01100  ; alpha CH2-N of sulfonamide. Use q=0.45 for CRH-N, q=0.65 for O=N-C-CH-N.
 opls_485    1.00800  ; H of alpha CH2-N of sulfonamide
 opls_486   12.01100  ; beta CH3 of N-ethyl sulfonamide
 opls_487    1.00800  ; H of beta CH3 of N-ethyl sulfonamide
 opls_488   12.01100  ; benzene C attached to S of sulfonamide
 opls_490   12.01100  ; C(H2OS) ethyl ester
 opls_491   12.01100  ; C(HOS) i-pr ester
 opls_492   12.01100  ; C(OS) t-bu ester
 opls_493   32.06000  ; S in sulfone    
 opls_494   15.99940  ; O in sulfone
 opls_496   32.06000  ; sulfoxide - all atom
 opls_497   15.99940  ; sulfoxide - all atom
 opls_498   12.01100  ; CH3 all-atom C: sulfoxide
 opls_499   12.01100  ; CH2 all-atom C: sulfoxide
 opls_500   12.01100  ; CG in Trp
 opls_501   12.01100  ; CD C in Trp
 opls_502   12.01100  ; CE C in Trp
 opls_503   14.00670  ; NE in Trp
 opls_504    1.00800  ; H on NE in Trp
 opls_505   12.01100  ; CB in His
 opls_506   12.01100  ; CE1 in HID, HIE
 opls_507   12.01100  ; CD2 in HID, CG in HIE
 opls_508   12.01100  ; CG in HID, CD2 in HIE
 opls_509   12.01100  ; CE1 in HIP
 opls_510   12.01100  ; CG, CD2 in HIP
 opls_511   14.00670  ; NE in HID, ND in HIE
 opls_512   14.00670  ; N in HIP
 opls_513    1.00800  ; H on N in HIP
 opls_514   12.01100  ; CD1 in TRP
 opls_515   12.01100  ; all-atom C: CH, isopropyl benzene
 opls_516   12.01100  ; all-atom C: C,  t-butyl benzene
 opls_517   12.01100  ; vinyl ether HCOR         
 opls_518   12.01100  ; vinyl ether RCOR         
 opls_520   14.00670  ; N   in pyridine 6-31G*
 opls_521   12.01100  ; C1  in pyridine CHELPG
 opls_522   12.01100  ; C2  in pyridine charges
 opls_523   12.01100  ; C3  in pyridine for
 opls_524    1.00800  ; H1  in pyridine 520-619
 opls_525    1.00800  ; H2  in pyridine
 opls_526    1.00800  ; H3  in pyridine
 opls_527   14.00670  ; N   in pyrazine
 opls_528   12.01100  ; C   in pyrazine
 opls_529    1.00800  ; H   in pyrazine
 opls_530   14.00670  ; N   in pyrimidine
 opls_531   12.01100  ; C2  in pyrimidine
 opls_532   12.01100  ; C4  in pyrimidine
 opls_533   12.01100  ; C5  in pyrimidine
 opls_534    1.00800  ; H2  in pyrimidine
 opls_535    1.00800  ; H4  in pyrimidine
 opls_536    1.00800  ; H5  in pyrimidine
 opls_537   14.00670  ; N   in pyridazine
 opls_538   12.01100  ; C3  in pyridazine
 opls_539   12.01100  ; C4  in pyridazine
 opls_540    1.00800  ; H3  in pyridazine
 opls_541    1.00800  ; H4  in pyridazine
 opls_542   14.00670  ; N   in pyrrole
 opls_543   12.01100  ; C2  in pyrrole
 opls_544   12.01100  ; C3  in pyrrole
 opls_545    1.00800  ; H1  in pyrrole
 opls_546    1.00800  ; H2  in pyrrole
 opls_547    1.00800  ; H3  in pyrrole
 opls_548   14.00670  ; N1  in pyrazole
 opls_549   14.00670  ; N2  in pyrazole
 opls_550   12.01100  ; C3  in pyrazole
 opls_551   12.01100  ; C4  in pyrazole
 opls_552   12.01100  ; C5  in pyrazole
 opls_553    1.00800  ; H1  in pyrazole
 opls_554    1.00800  ; H3  in pyrazole
 opls_555    1.00800  ; H4  in pyrazole
 opls_556    1.00800  ; H5  in pyrazole
 opls_557   14.00670  ; N1  in imidazole
 opls_558   12.01100  ; C2  in imidazole
 opls_559   14.00670  ; N3  in imidazole
 opls_560   12.01100  ; C4  in imidazole
 opls_561   12.01100  ; C5  in imidazole
 opls_562    1.00800  ; H1  in imidazole
 opls_563    1.00800  ; H2  in imidazole
 opls_564    1.00800  ; H4  in imidazole
 opls_565    1.00800  ; H5  in imidazole
 opls_566   15.99940  ; O   in furan
 opls_567   12.01100  ; C2  in furan
 opls_568   12.01100  ; C3  in furan
 opls_569    1.00800  ; H2  in furan
 opls_570    1.00800  ; H3  in furan
 opls_571   15.99940  ; O   in oxazole
 opls_572   12.01100  ; C2  in oxazole
 opls_573   14.00670  ; N   in oxazole
 opls_574   12.01100  ; C4  in oxazole
 opls_575   12.01100  ; C5  in oxazole
 opls_576    1.00800  ; H2  in oxazole
 opls_577    1.00800  ; H4  in oxazole
 opls_578    1.00800  ; H5  in oxazole
 opls_579   15.99940  ; O   in isoxazole
 opls_580   14.00670  ; N   in isoxazole
 opls_581   12.01100  ; C3  in isoxazole
 opls_582   12.01100  ; C4  in isoxazole
 opls_583   12.01100  ; C5  in isoxazole
 opls_584    1.00800  ; H3  in isoxazole
 opls_585    1.00800  ; H4  in isoxazole
 opls_586    1.00800  ; H5  in isoxazole
 opls_587   14.00670  ; N1  in indole
 opls_588   12.01100  ; C2  in indole
 opls_589   12.01100  ; C3  in indole
 opls_590   12.01100  ; C4  in indole
 opls_591   12.01100  ; C5  in indole
 opls_592   12.01100  ; C6  in indole
 opls_593   12.01100  ; C7  in indole
 opls_594   12.01100  ; C8  in indole
 opls_595   12.01100  ; C9  in indole
 opls_596    1.00800  ; H1  in indole
 opls_597    1.00800  ; H2  in indole
 opls_598    1.00800  ; H3  in indole
 opls_599    1.00800  ; H4  in indole
 opls_600    1.00800  ; H5  in indole
 opls_601    1.00800  ; H6  in indole
 opls_602    1.00800  ; H7  in indole
 opls_603   14.00670  ; N1  in quinoline
 opls_604   12.01100  ; C2  in quinoline
 opls_605   12.01100  ; C3  in quinoline
 opls_606   12.01100  ; C4  in quinoline
 opls_607   12.01100  ; C5  in quinoline
 opls_608   12.01100  ; C6  in quinoline
 opls_609   12.01100  ; C7  in quinoline
 opls_610   12.01100  ; C8  in quinoline
 opls_611   12.01100  ; C9  in quinoline
 opls_612   12.01100  ; C10 in quinoline
 opls_613    1.00800  ; H2  in quinoline
 opls_614    1.00800  ; H3  in quinoline
 opls_615    1.00800  ; H4  in quinoline
 opls_616    1.00800  ; H5  in quinoline
 opls_617    1.00800  ; H6  in quinoline
 opls_618    1.00800  ; H7  in quinoline
 opls_619    1.00800  ; H8  in quinoline
 opls_620   14.00670  ; N1  in purine 
 opls_621   12.01100  ; C2  in purine   
 opls_622   14.00670  ; N3  in purine   
 opls_623   12.01100  ; C4  in purine
 opls_624   12.01100  ; C5  in purine
 opls_625   12.01100  ; C6  in purine
 opls_626   14.00670  ; N7  in purine
 opls_627   12.01100  ; C8  in purine
 opls_628   14.00670  ; N9  in purine
 opls_629    1.00800  ; H2  in purine
 opls_630    1.00800  ; H6  in purine
 opls_631    1.00800  ; H8  in purine
 opls_632    1.00800  ; H9  in purine
 opls_633   32.06000  ; S   in thiazole
 opls_634   12.01100  ; C2  in thiazole
 opls_635   14.00670  ; N   in thiazole
 opls_636   12.01100  ; C4  in thiazole
 opls_637   12.01100  ; C5  in thiazole
 opls_638    1.00800  ; H2  in thiazole
 opls_639    1.00800  ; H4  in thiazole
 opls_640    1.00800  ; H5  in thiazole
 opls_641   14.00670  ; N   in 1,3,5-triazine
 opls_642   12.01100  ; C   in 1,3,5-triazine
 opls_643    1.00800  ; H   in 1,3,5-triazine
 opls_644   12.01100  ; C5  in serotonin
 opls_645   12.01100  ; C on C3 in serotonin
 opls_646   14.00670  ; N1,N10   in 1,10-phenanthroline
 opls_647   12.01100  ; C2,C9  in 1,10-phenanthroline
 opls_648   12.01100  ; C3,C8  in 1,10-phenanthroline
 opls_649   12.01100  ; C4,C7  in 1,10-phenanthroline
 opls_650   12.01100  ; C12,C14 in 1,10-phenanthroline
 opls_651   12.01100  ; C11,C13 in 1,10-phenanthroline
 opls_652   12.01100  ; C5  in 1,10-phenanthroline
 opls_653    1.00800  ; H2,H9  in 1,10-phenanthroline
 opls_654    1.00800  ; H3,H8  in 1,10-phenanthroline
 opls_655    1.00800  ; H4,H7  in 1,10-phenanthroline
 opls_656    1.00800  ; H5,H6  in 1,10-phenanthroline
 opls_670   12.01100  ; CH3, 2-methyl pyridine
 opls_671   12.01100  ; CH2, 2-ethyl pyridine
 opls_672   12.01100  ; CH3, 3-methyl pyridazine
 opls_673   12.01100  ; CH2, 3-ethyl pyridazine
 opls_674   12.01100  ; CH3, 4-methyl pyrimidine
 opls_675   12.01100  ; CH2, 4-ethyl pyrimidine
 opls_676   12.01100  ; CH3, 2-methyl pyrazine
 opls_677   12.01100  ; CH2, 2-ethyl pyrazine
 opls_678   12.01100  ; CH3, 2-methyl pyrrole
 opls_679   12.01100  ; CH2, 2-ethyl pyrrole
 opls_680   12.01100  ; CH3, 2-methyl furan
 opls_681   12.01100  ; CH2, 2-ethyl furan
 opls_697    0.00000  ; Ac+3 Actinide params -
 opls_698    0.00000  ; Th+4
 opls_699    0.00000  ; Am+3 F. van Veggel
 opls_700   12.01100  ; C+  in t-butyl+ B3LYP/6-31G*
 opls_701   12.01100  ; C   in t-butyl+   charges
 opls_702    1.00800  ; H   in t-butyl+
 opls_703    0.00000  ; La+3
 opls_704    0.00000  ; Nd+3 Lanthanide params -
 opls_705    0.00000  ; Eu+3 F. van Veggel, Chem Eur J 5, 90 (1999).
 opls_706    0.00000  ; Gd+3             
 opls_707    0.00000  ; Yb+3 see also JPC-A 104, 7659 (2000)
 opls_708   12.01100  ; C  in Cl..CH3..Cl- TS
 opls_709   35.45300  ; Cl charges: JACS 117,2024 (95)
 opls_710    1.00800  ; H  in Cl..CH3..Cl- TS
 opls_711   12.01100  ; CH2    C: cyclopropane
 opls_712   12.01100  ; CHR    C: cyclopropane
 opls_713   12.01100  ; CR2    C: cyclopropane
 opls_714   12.01100  ; C in C5H5- cyclopentadienyl anion
 opls_715    1.00800  ; H in C5H5- cyclopentadienyl anion
 opls_716   12.01100  ; C in C5H5  cyclopentadienyl radical
 opls_717    1.00800  ; H in C5H5  cyclopentadienyl radical
 opls_718   12.01100  ; C(F)  fluorobenzene
 opls_719   18.99840  ; F     fluorobenzene
 opls_720   12.01100  ; C(F)  hexafluorobenzene
 opls_721   18.99840  ; F     hexafluorobenzene
 opls_722   79.90400  ; Br    alkyl bromide (UA, but probably ok for AA)
 opls_724   12.01100  ; C(CF3) trifluoromethylbenzene
 opls_725   12.01100  ; CF3   trifluoromethylbenzene
 opls_726   18.99840  ; F     trifluoromethylbenzene
 opls_727   12.01100  ; C(F)  difluorobenzenes
 opls_728   18.99840  ; F     difluorobenzenes
 opls_729   12.01100  ; C(Br) bromobenzene
 opls_730   79.90400  ; Br    bromobenzene
 opls_731   12.01100  ; C(I)  iodobenzene - tentative
 opls_732  126.90450  ; I     iodobenzene - tentative
 opls_733   12.01100  ; all-atom C: CH, cyclopropyl benzene
 opls_734   32.06000  ; all-atom S: thiophenol (HS is #204)
 opls_735   12.01100  ; C(S)  thiophenol
 opls_736   12.01100  ; CG of Benzamidine
 opls_737   12.01100  ; CD of Benzamidine
 opls_738   12.01100  ; CE of Benzamidine
 opls_739   12.01100  ; CZ of Benzamidine
 opls_740    1.00800  ; HD of Benzamidine
 opls_741    1.00800  ; HE of Benzamidine
 opls_742   12.01100  ; C+ of Benzamidine
 opls_743   14.00670  ; N-H2 of Benzamidine
 opls_744    1.00800  ; H1-N of Benzamidine
 opls_745    1.00800  ; H2-N of Benzamidine
 opls_746    1.00800  ; H-CG of Benzamidine
 opls_747   12.01100  ; CH3 in neutral MeGDN
 opls_748   12.01100  ; CD of neutral ARG
 opls_749   14.00670  ; NE of neutral ARG
 opls_750   14.00670  ; N1 of neutral ARG (HN=CZ)
 opls_751   14.00670  ; N2 of neutral ARG (H2N-CZ)
 opls_752   12.01100  ; CZ of neutral ARG
 opls_753   14.00670  ; N IN RCN  nitriles
 opls_754   12.01100  ; C IN RCN  nitriles
 opls_755   12.01100  ; C of CH3 in  CH3CN
 opls_756   12.01100  ; C of CH2 in RCH2CN
 opls_757   12.01100  ; C of CH  in R2CHCN
 opls_758   12.01100  ; C of C   in R3CCN
 opls_759    1.00800  ; HC-CT-CN alpha-H in nitriles
 opls_760   14.00670  ; N in nitro R-NO2
 opls_761   15.99940  ; O in nitro R-NO2
 opls_762   12.01100  ; CT-NO2 nitromethane
 opls_763    1.00800  ; HC-CT-NO2 alpha-H in nitroalkanes
 opls_764   12.01100  ; CT-NO2 nitroethane
 opls_765   12.01100  ; CT-NO2 2-nitropropane
 opls_766   12.01100  ; CT-NO2 2-methyl-2-nitropropane
 opls_767   14.00670  ; N in nitro Ar-NO2
 opls_768   12.01100  ; C(NO2) nitrobenzene
 opls_771   15.99940  ; propylene carbonate O (Luciennes param.)
 opls_772   12.01100  ; propylene carbonate C=O   
 opls_773   15.99940  ; propylene carbonate OS    
 opls_774   12.01100  ; propylene carbonate C in CH2
 opls_775   12.01100  ; propylene carbonate C in CH
 opls_776   12.01100  ; propylene carbonate C in CH3
 opls_777    1.00800  ; propylene carbonate H in CH2
 opls_778    1.00800  ; propylene carbonate H in CH
 opls_779    1.00800  ; propylene carbonate H in CH3
 opls_781   30.97376  ; phosphonium R4P+
 opls_782   12.01100  ; CH3PR3+ 6-31G* CHELPG
 opls_783   12.01100  ; RCH2PR3+
 opls_784    1.00800  ; H in CH3PR3+
 opls_785   30.97376  ; P in PF6-
 opls_786   18.99840  ; F in PF6-
 opls_787   14.00670  ; N in NO3- 
 opls_788   15.99940  ; O in NO3- 
 opls_795   15.99940  ; O TIP4F Water  
 opls_796    1.00800  ; H TIP4F Water  
 opls_797    0.00000  ; M TIP4F Water  
 opls_900   14.00670  ; N primary   amines
 opls_901   14.00670  ; N secondary amines, aziridine N1 
 opls_902   14.00670  ; N tertiary  amines
 opls_903   12.01100  ; CH3(N) primary   aliphatic amines, H(C) is #911
 opls_904   12.01100  ; CH3(N) secondary aliphatic amines, H(C) is #911
 opls_905   12.01100  ; CH3(N) tertiary  aliphatic amines, H(C) is #911
 opls_906   12.01100  ; CH2(N) primary   aliphatic amines, H(C) is #911
 opls_906B  12.01100  ; CA in GLY-NH2 N-terminus
 opls_907   12.01100  ; CH2(N) secondary aliphatic amines, aziridine  C2,C3H
 opls_908   12.01100  ; CH2(N) tertiary  aliphatic amines, H(C) is #911
 opls_909    1.00800  ; H(N)   primary   amines
 opls_910    1.00800  ; H(N)   secondary amines
 opls_911    1.00800  ; H(C) for C bonded to N in amines, diamines (aziridine H2,H3)
 opls_912   12.01100  ; CH     primary isopropyl amine
 opls_912B  12.01100  ; CA in NH2 N-terminus. All AA except GLY, PRO
 opls_913   12.01100  ; C      primary t-butyl amine
 opls_914   12.01100  ; CH     secondary isopropyl amine
 opls_915   12.01100  ; CH     tertiary  isopropyl amine
 opls_916   12.01100  ; C(NH2) aniline
 opls_917   12.01100  ; C(NH2) N-methylaniline
 opls_918   12.01100  ; C(NH2) N,N-dimethylaniline
 opls_925   12.01100  ; alkyne RC%CH terminal C   acetylene
 opls_926    1.00800  ; alkyne RC%CH terminal H
 opls_927   12.01100  ; alkyne RC%CH C2 R-with 2 or 3 H
 opls_928   12.01100  ; alkyne RC%CH C2 R-with 1 H
 opls_929   12.01100  ; alkyne RC%CH C2 R-with no H or R=Phenyl
 opls_930    1.00800  ; alkyne RC%CH H on C3 (for C3 use #135-#139)
 opls_931   12.01100  ; alkyne RC%CR
 opls_940   14.00670  ; N (R3NH+)
 opls_941    1.00800  ; H (R3NH+)
 opls_942   12.01100  ; C in  CH3NHR2+
 opls_943   12.01100  ; C in  RCH2NHR2+
 opls_944   12.01100  ; C in  R2CHNHR2+
 opls_945   12.01100  ; C in  R3CNHR2+
 opls_950    1.00800  ; glycine zwit. 6-31G* CHELPG charges
 opls_951   12.01100  ; glycine zwit. 6-31G* CHELPG charges
 opls_952   12.01100  ; glycine zwit. 6-31G* CHELPG charges
 opls_953   14.00670  ; glycine zwit. 6-31G* CHELPG charges
 opls_954   15.99940  ; glycine zwit. 6-31G* CHELPG charges
 opls_955    1.00800  ; glycine zwit. 6-31G* CHELPG charges
 opls_956   18.99840  ; F  in monoalkyl fluorides (tentative)
 opls_957   12.01100  ; RCH2F in monoalkyl fluorides (tentative)
 opls_958    1.00800  ; H in RCHF in monoalkyl fluorides (tentative)
 opls_959   12.01100  ; R2CHF in monoalkyl fluorides (tentative)
 opls_960   12.01100  ; R3CF in monoalkyl fluorides (tentative)
 opls_961   12.01100  ; CF3 perfluoroalkanes
 opls_962   12.01100  ; CF2 perfluoroalkanes
 opls_963   12.01100  ; CF perfluoroalkanes
 opls_964   12.01100  ; CF4
 opls_965   18.99840  ; F: perfluoroalkanes
 MNH3        0.0      ; Dummy mass in rigid tetraedrical NH3 group
 MNH2        0.0      ; Dummy mass in rigid umbrella-shaped NH2 group
 MCH3A       0.0      ; Dummy mass in rigid tetraedrical CH3 group
 MCH3B       0.0      ; Dummy mass in rigid tetraedrical CH3 group
 MW          0.0      ; Dummy mass in rigid tyrosine rings
 DUM         0.0      ; Dummy mass in TIP4P etc.
; These ion atomtypes are NOT part of OPLS, but since they are
; needed for some proteins we have added them.
 Cu2+       63.546    ; Copper. See Inorg. Chem. 40, 5223 (2001).
 Fe2+       55.847    ; Iron 
 Zn2+       65.370    ; Zinc 
 Ar         39.948    ; Argon
; Added by DvdS 05/2005 copied from GROMACS force field.       
 SI         28.080    ; Silicium in Glass etc.
"""
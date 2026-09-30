# shifts in units of 1/L (alpha); (a,b,c) triples actually used (coincident ones perturbed +-2%)
D=0.02
TRIPLES=[  # (label, nominal, actual)
 ('(1,1,1)',(1,1,1),(1-D,1+D,1.0)),
 ('(2,1,0.5)',(2,1,0.5),(2.0,1.0,0.5)),
 ('(0.5,2,1)',(0.5,2,1),(0.5,2.0,1.0)),
 ('(3,3,2)',(3,3,2),(3*(1-D),3*(1+D),2.0)),
 ('(1,4,0.5)',(1,4,0.5),(1.0,4.0,0.5)),
 ('(0.3,0.6,0.4)',(0.3,0.6,0.4),(0.3,0.6,0.4)),
 ('(2,2,4)',(2,2,4),(2*(1-D),2*(1+D),4.0)),
]
PAIRS_XXB=[(0.5,0.5),(1.0,1.0),(2.0,2.0),(1.0,2.0),(0.5,2.0),(3.0,3.0),(0.3,0.6),(1.0,4.0)]
PAIRS_XX=[(1.0,2.0),(0.5,2.0),(1.0,4.0),(0.3,0.6),(0.5,1.0),(2.0,4.0)]
def alpha_list():
    s=set()
    for _,_,t in TRIPLES: s.update(round(x,10) for x in t)
    for p in PAIRS_XXB+PAIRS_XX: s.update(round(x,10) for x in p)
    return sorted(s)

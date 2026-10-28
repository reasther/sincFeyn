import re, sys, os
def load(g):
    src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'k', f'{g}_k.py')).read()
    src=re.sub(r'iV = 1/\((.*)\)', r'vv = (\1)', src)
    src=src.replace('ivp = 1\n','nn = 0\n').replace('ivp = ivp*iV','nn += 1; ivp = PW(vv, nn)')
    src=src.replace('def pt(x):','def pt(x, PW):')
    src=src.replace('*pt(s)','*pt(s, PW)').replace('def f(sf):','def f(sf, PW):')
    ns={}; exec(src,ns); return ns

# This example assumes we have a mesh object selected

import bpy
import bmesh
import numpy as np

def v2pDistance(v, neig_v):
    vp_co=v.co
    v.normal_update()
    n = v.normal
    neig_v_co = neig_v.co
    dist = (neig_v_co-vp_co)
    t = dist.length
    h = n.dot(dist)
    return t,h

#    return exp(-(t*t)/2*p0*p0),exp(-(h*h)/2*p1*p1)
    

def neighborhood(v, vcallback):
    elinks = v.link_edges
    th = []
    for e in elinks:
        neig_v = e.other_vert(v)
        th.append(vcallback(v,neig_v))
    return v,np.array(th)

def averageEdgeLenght(edges):
    sum = 0
    for e in edges:
        sum += e.calc_length()
    return sum/len(edges)
    
def printVec(vec):
    for v in vec:
        print(v.shape)
#        print(v)

def gaussion(x,sigma):
    return np.exp(-0.5*(x/sigma)**2)
    
def Iter(bm):
    d = []
    for v in bm.verts:
        d.append(neighborhood(v, v2pDistance))
    mu0 = averageEdgeLenght(bm.edges)
    mu0 = mu0*1.4
#    return d
#    print("Iter")
    if mu0 != 0:
        for v,info in d:
            h = info[:,1]
            mu1 = h.mean()
            if mu1 != 0:
                t = info[:,0]
                wc = gaussion(t,mu0)
                ws = gaussion(h,mu1)
                if wc.sum() > 10E-5 and ws.sum() > 10E-5:
                    print(t.shape,h.shape,mu0,mu1)
                    normalizer = wc*ws
                    sum = normalizer * h
                    v.co = v.co + v.normal*(sum.sum()/normalizer.sum())
                
def normal_distribution(x, mean, sigma):
    return np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi)* sigma)

def getBM():
    # Get the active mesh
    me = bpy.context.object.data
    # Get a BMesh representation
    bm = bmesh.new()   # create an empty BMesh
    bm.from_mesh(me)   # fill it in from a Mesh
    return bm

bm = getBM()

for x in range(5):
    Iter(bm)


# Finish up, write the bmesh back to the mesh
bm.to_mesh(bpy.context.object.data)
bm.free()  # free and prevent further access
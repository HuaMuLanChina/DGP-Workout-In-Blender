import bpy
import bmesh
import numpy as np

def getBM():
    # Get the active mesh
    me = bpy.context.object.data
    # Get a BMesh representation
    bm = bmesh.new()   # create an empty BMesh
    bm.from_mesh(me)   # fill it in from a Mesh
    return bm

def localaxis(v0, v1, v2):
    vec0 = v1.co - v0.co
    vec1 = v2.co - v0.co
    x = vec0.normalized()
    n = (x.cross(vec1)).normalized()
    y = n.cross(x)
    return [vec0.length,vec1.dot(x),vec1.dot(y)]

def projected2D(fs, vs):
    fs_num = len(fs)
    d = np.zeros((fs_num,3,2))
    for f in fs:
        d0,d1,d2 = localaxis(f.verts[0],f.verts[1],f.verts[2])
        f_i = f.index
        #d[f_i][0][0] = 0
        #d[f_i][0][1] = 0
        d[f_i][1][0] = d0
        #d[f_i][1][1] = 0
        d[f_i][2][0] = d1
        d[f_i][2][1] = d2
    return d

def buildM(bm):
    fs_num = len(bm.faces)
    vs_num = len(bm.verts)
    d = projected2D(bm.faces, bm.verts)
    M = np.zeros((fs_num, vs_num, 2)) #((fs_num, vs_num, real&img))
    for f in bm.faces:
        f_i = f.index
        area = np.sqrt(f.calc_area())
        M[f_i,[f.verts[0].index, f.verts[1].index, f.verts[2].index],:] = (d[f_i,[2,0,1],:] - d[f_i,[1,2,0],:]) / area
    free,fix = M[:,:-2,:],M[:,-2:,:]
    real,img = free[...,0],free[...,1]
    fixReal,fixImg = fix[...,0],fix[...,1]
    return real,img,fixReal,fixImg

def composeC(r,i):
    t0 = np.concatenate([r,-i], axis=1)
    t1 = np.concatenate([i,r], axis=1)
    t = np.concatenate([t0,t1])
    return t
bm = getBM()
r,i,fr,fi = buildM(bm)
A = composeC(r,i)
B = composeC(fr,fi)
C = -B@np.array([[0],[1],[0],[0]])
At = A.transpose()
x = np.linalg.solve(At@A,At@C)

#num_verts = len(bm.verts)
#for i,v in enumerate(bm.verts):
#    v.co.z = 0
#    if i < num_verts -2:
#        v.co.x = x[i]
#        v.co.y = x[i+num_verts-2]
#    elif i == num_verts - 2:
#        v.co.x = 0
#        v.co.y = 0
#    else:
#        v.co.x = 1
#        v.co.y = 0
#normallize to [0.1]
uv = x.reshape(-1, x.size//2).transpose()
fixuv = np.array([[0,0],[1,0]])
uv_all = np.vstack([uv,fixuv])
orig = uv_all.min(axis=0)
trans = np.array([0,0]) - orig
uv_all = uv_all + trans
uv_all = uv_all/uv_all.max()
#set to uv channel
me = bpy.context.object.data
uv_layer = me.uv_layers.active.data
for poly in me.polygons:
    for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
        uv_layer[loop_index].uv.x = uv_all[me.loops[loop_index].vertex_index,0]
        uv_layer[loop_index].uv.y = uv_all[me.loops[loop_index].vertex_index,1]

bm.free()
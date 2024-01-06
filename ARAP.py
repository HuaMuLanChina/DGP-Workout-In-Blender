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

def buidlFV(fs,vs,det):
    fs_num = len(fs)
    vs_num = len(vs)
    FV = np.zeros((fs_num*4, vs_num*2))
    for f in fs:
        fi = f.index
        vii,vji,vki = f.vertices
        FV[fi*4  ,[vii,vji,vki]] = det[fi,0]
        FV[fi*4+1,[vii,vji,vki]] = det[fi,1]
        FV[fi*4+2,[vii+vs_num,vji+vs_num,vki+vs_num]] = det[fi,0]
        FV[fi*4+3,[vii+vs_num,vji+vs_num,vki+vs_num]] = det[fi,1]
    return FV

def localaxis(v0, v1, v2):
    vec0 = v1.co - v0.co
    vec1 = v2.co - v0.co
    x = vec0.normalized()
    n = (x.cross(vec1)).normalized()
    y = n.cross(x)
    return [vec0.length,vec1.dot(x),vec1.dot(y)]

def calDet(fs, vs):
    fs_num = len(fs)
    et = np.zeros((fs_num,2,3))
    for f in fs:
        d = np.zeros((3,2))
        vii,vji,vki = f.vertices
        d0,d1,d2 = localaxis(vs[vii],vs[vji],vs[vki])
        #d[f_i][0][0] = 0
        #d[f_i][0][1] = 0
        d[1,0] = d0
        #d[f_i][1][1] = 0
        d[2,0] = d1
        d[2,1] = d2
        vi,vj,vk = d[0],d[1],d[2]
        a = 1/(2*f.area)
        et[f.index] = a * np.array([[vj[1]-vk[1], vk[1]-vi[1], vi[1]-vj[1]],[vk[0]-vj[0],vi[0]-vk[0],vj[0]-vi[0]]])
    return et

def get2Duv(me):
    num_vs = len(me.vertices)
    uv_layer = me.uv_layers.active.data
    uv = np.zeros((num_vs,2))
    for poly in me.polygons:
        for loop_index in poly.loop_indices:
            vi = me.loops[loop_index].vertex_index
            uv[vi,0] = uv_layer[loop_index].uv.x
            uv[vi,1] = uv_layer[loop_index].uv.y
    return uv
            
def set2Duv(me, uv):
    #Rescale to range[0,1] begin
    orig = uv.min(axis=0)
    trans = np.array([0,0]) - orig
    uv = uv + trans
    uv = uv/uv.max()
    #Rescale to range[0,1] end
    uv_layer = me.uv_layers.active.data
    for poly in me.polygons:
        for loop_index in poly.loop_indices:
            vi = me.loops[loop_index].vertex_index
            uv_layer[loop_index].uv.x = uv[vi,0]
            uv_layer[loop_index].uv.y = uv[vi,1]

def localOpt(fs,vs,det,uv):
    fs_num = len(fs)
    Lt = np.zeros((fs_num,2,2))
    E = 0
    for f in fs:
        fi = f.index
        vii,vji,vki = f.vertices
        et = det[fi]
        fuv = np.array([uv[vii],uv[vji],uv[vki]])
        #print(vi.shape,et.shape,fuv.shape)
        jac = (et@fuv).transpose()
        U, S, Vh = np.linalg.svd(jac)
        if np.linalg.det(U@Vh) > 0:
            Lt[fi] = U@Vh
        else:
            print("flatten")
            minus = np.identity(2)
            minus[1,1] = -1
            Lt[fi] = U@(Vh@minus)
        E += np.sum((jac-Lt[fi])**2)
    return Lt,E

me = bpy.context.object.data
uv = get2Duv(me)
polys = me.polygons
verts = me.vertices
det = calDet(polys,verts)
FV = buidlFV(polys,verts,det)
FVt = FV.transpose()
Lt,E = localOpt(polys,verts,det,uv)
for i in range(5):
    print(f"Energe:{E}")
    newUV = np.linalg.solve(FVt@FV,FVt@Lt.flatten()).reshape(2,-1).transpose()
    Lt,E = localOpt(polys,verts,det,newUV)

set2Duv(me, newUV)
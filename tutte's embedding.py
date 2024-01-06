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

#find open boundary , multi boundary find the one with max number edges   BENDG
def is_include(e, exclude):
    for ex in exclude:
        if e in ex:
            return True
    return False

def findOpenEdge(edges, exclude):
    for e in edges:
        if e.is_boundary and not is_include(e, exclude):
            return e
    return None

def nextBoundaryEdge(e,preE):
    stv,nxv = e.verts
    for nxe in nxv.link_edges:
        if nxe.is_boundary and nxe != e and nxe != preE:
            return nxe
    for nxe in stv.link_edges:
        if nxe.is_boundary and nxe != e and nxe != preE:
            return nxe

def boundaryloop(e):
    print(f"start:{e}")
    boundary = [e]
    nxe = nextBoundaryEdge(e, None)
    preE = e
    while nxe != e:
        print("---")
        print(nxe,type(nxe))
        boundary.append(nxe)
        nxe = nextBoundaryEdge(boundary[-1], preE)
        preE = boundary[-1]
    return boundary

def findBoundary(bm):
    e = findOpenEdge(bm.edges, [])
    bds = []
    while e != None:
        bds.append(boundaryloop(e))
        e = findOpenEdge(bm.edges, bds)
    return bds
#find open boundary , multi boundary find the one with max number edges END

def boundaryVerts(boundaryEdges):
    stv,nxv = boundaryEdges[0].verts
    verts = [stv]
    for e in boundaryEdges[0:-1]:
        other = e.other_vert(stv)
        verts.append(other)
        stv = other
    return verts

def bdE2bdV(bdE):
    bdv = []
    for bde in bdE:
        bdv.append(boundaryVerts(bde))
    return bdv

#tutte put into a circle boundary in order
def circleBoundary(bdV,r):
    delta = 2 * np.pi / len(bdV)
    for i,v in enumerate(bdV):
        v.co.z = 0
        v.co.x = r*np.cos(delta*i)
        v.co.y = r*np.sin(delta*i)

def buildLaplacianB(bdV, verts):
    n = len(verts)
    Lx = np.zeros((n,n))
    Bx = np.zeros(n)
    Ly = np.zeros((n,n))
    By = np.zeros(n)
    #use bmesh vert index for indexing
    for v in verts:
        i = v.index
        if v in bdV:
            Lx[i][i] = 1.
            Bx[i] = v.co.x
            Ly[i][i] = 1.
            By[i] = v.co.y
        else:
            Lx[i][i] = -1.
            Bx[i] = 0.
            Ly[i][i] = -1.
            By[i] = 0.
            neigbs = [e.other_vert(v) for e in v.link_edges]
            idegree = 1/len(neigbs)
            for neigb in neigbs:
                ni = neigb.index
                Lx[i][ni] = Ly[i][ni] = idegree           
    return Lx,Bx,Ly,By

def setVerts(verts,x,y):
    for v in verts:
        i = v.index
        v.co.x = x[i]
        v.co.y = y[i]
        v.co.z = 0
    
bm = getBM()
bdE = findBoundary(bm)
bdV = bdE2bdV(bdE)
#multi boundary find the one with max number edges BENDG
max = bdV[0]
if len(bdV) > 1:
    for bd in bdV[1:]:
        if len(max) < len(bd):
            max = bd
#multi boundary find the one with max number edges END
circleBoundary(max,0.5)
Lx,Bx,Ly,By = buildLaplacianB(max,bm.verts)
x = np.linalg.solve(Lx,Bx)
y = np.linalg.solve(Ly,By)
x += 0.5
y += 0.5
#setVertsU(bm.verts,x,y)
#set uv coords  to mesh uv channel
me = bpy.context.object.data
uv_layer = me.uv_layers.active.data
for poly in me.polygons:
    #print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
    # range is used here to show how the polygons reference loops,
    # for convenience 'poly.loop_indices' can be used instead.
    for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
        uv_layer[loop_index].uv.x = x[me.loops[loop_index].vertex_index]
        uv_layer[loop_index].uv.y = y[me.loops[loop_index].vertex_index]
        #print("    Vertex: %d" % me.loops[loop_index].vertex_index)
        #print("    UV: %r" % uv_layer[loop_index].uv)
#bm.to_mesh(bpy.context.object.data)
bm.free()
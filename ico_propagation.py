# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 17:03:51 2017

@author: tigan_5ytncvu
"""
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np


t = (1.0 + (5.0)**(0.5) )/ 2.0

icosahedron_vertices = [[-1,  t,  0],
                            [ 1,  t,  0],
                            [-1, -t,  0],
                            [1, -t,  0],
                            [ 0, -1,  t],
                            [0,  1,  t],
                            [ 0, -1, -t],
                            [ 0,  1, -t],
                            [ t,  0, -1],
                            [t,  0,  1],
                            [-t,  0, -1],
                            [-t,  0,  1]]

faces = np.array(
            [[0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1]])



class Sphere:
    
    def __init__(self, vertices = icosahedron_vertices, faces = faces, recursion_level = 4 ):
        
        self.icosahedron_vertices = vertices
        self.faces = faces
        self.recursion_level = recursion_level
        self.ch = []


    
    def normalize_v3(self, arr):
        ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    
        lens = arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 
        lens = np.sqrt(lens)
        arr[:,0] /= lens
        arr[:,1] /= lens
        arr[:,2] /= lens        
        return arr
        
    def normalize(self, arr):
        ''' Normalize components '''
     
        lens = arr[0]**2 + arr[1]**2 + arr[2]**2 
        lens = np.sqrt(lens)
        arr[0] /= lens
        arr[1] /= lens
        arr[2] /= lens          
        return arr
    
    
    def splitter(self, vertices, faces):
    
        new_faces = []
        vert_list = []
        face_counter = 0
        
        for face in faces: 
            face_counter += 1
            #print("facecount =", face_counter)
            v0 = face[0]
            v1 = face[1]
            v2 = face[2]
            
            a = ( v0+v2 ) * 0.5
            b = ( v0+v1 ) * 0.5
            c = ( v1+v2 ) * 0.5  
            
            self.normalize( a )
            self.normalize( b )
            self.normalize( c )
    
            vert_list.append(v0)   
            vert_list.append(v1)
            vert_list.append(v2)
            vert_list.append(a)
            vert_list.append(b)
            vert_list.append(c)
            #if v0marker + bmarker + amarker == 0:
            new_faces.append([v0, b, a])
            #if bmarker + v1marker + cmarker == 0:
            new_faces.append([b,v1,c])
            #if amarker + bmarker + cmarker == 0
            new_faces.append([a,b,c])
            new_faces.append([a,c,v2])
    
        new_faces = np.asarray(new_faces)
        vert_list = np.asarray(vert_list)
    
        return  vert_list, new_faces
    
    
    def create_unit_sphere_vert(self):
        vertex_array, index_array = self.icosahedron_vertices, self.faces
        for i in range( self.recursion_level - 1 ):
            vertex_array, index_array = self.splitter(vertex_array, index_array)
        return vertex_array, index_array
    
    
    def sph_polar_convert(self, xyz):
        """converts xyz array into r,theta,phi array"""
        ptsnew = np.zeros(xyz.shape)
        xy = xyz[:,0]**2 + xyz[:,1]**2
        ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
        ptsnew[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
        #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
        ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
        return ptsnew  
    
    def unique_rows(self, a):
        a = np.ascontiguousarray(a)
        unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
        return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    
    def unique(self, a):
        order = np.lexsort(a.T)
        a = a[order]
        diff = np.diff(a, axis=0)
        ui = np.ones(len(a), 'bool')
        ui[1:] = (diff != 0).any(axis=1) 
        return a[ui]
        
    def merge(self, a, b):
        max_offset = len(b)  # can't overlap with greater size than len(b)
        for i in reversed(range(max_offset+1)):
            # checks for equivalence of decreasing sized slices
            if a[-i:] == b[:i]:
                break
        return a + b[i:]
        
    def remove_duplicates(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    def remove_unique(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if (x in seen or seen_add(x))]
        
        
    def facet_ncheck(self,  facet_index, common_pt):#, f_list):
                
        nind = self.ch.neighbors[facet_index] #indices of the neighbouring facets
        facet1, facet2, facet3 = self.ch.simplices[nind[0]], self.ch.simplices[nind[1]], self.ch.simplices[nind[2]]  #facets 1,2 & 3
        adj_facet_list = []
        neigh_ind = []
        
        if (np.all(self.ch.points[facet1[0]] == common_pt) or 
                  np.all(self.ch.points[facet1[1]] == common_pt) or 
                         np.all(self.ch.points[facet1[2]] == common_pt)):
            adj_facet_list.append(facet1)
            neigh_ind.append(nind[0])
            
        if (np.all(self.ch.points[facet2[0]] == common_pt) or 
                  np.all(self.ch.points[facet2[1]] == common_pt) or 
                         np.all(self.ch.points[facet2[2]] == common_pt)):
            adj_facet_list.append(facet2)
            neigh_ind.append(nind[1])
            
        if (np.all(self.ch.points[facet3[0]] == common_pt) or 
                  np.all(self.ch.points[facet3[1]] == common_pt) or 
                         np.all(self.ch.points[facet3[2]] == common_pt)):
            adj_facet_list.append(facet3)
            neigh_ind.append(nind[2])

        return adj_facet_list, neigh_ind
        
        
            
    def find_pentagon(self):       
        
        for i in range(len(self.ch.simplices)):
            for point in self.ch.simplices[i]: 
                if np.all(self.ch.points[point] == self.icosahedron_vertices[0]): #means that point belongs to a pentaagon
                    s_facet = self.ch.simplices[i]
                    facet_nlist, neigh_ind = self.facet_ncheck( facet_index =i, common_pt =self.icosahedron_vertices[0])
                    while len(neigh_ind) < 6:
                        if len(neigh_ind) == 5:
                            facet_nlist.append(s_facet)
                            return facet_nlist, neigh_ind
                        neigh_facets, neigh_indices = self.facet_ncheck( facet_index =neigh_ind[-1], common_pt =self.icosahedron_vertices[0] )
                        neigh_ind = self.remove_duplicates(self.merge(list(neigh_ind),list(neigh_indices)))
                    return facet_nlist, neigh_ind
                    
    def find_hex(self, commonpt):       
        
        for i in range(len(self.ch.simplices)):
            for point in self.ch.simplices[i]: 
                if np.all(self.ch.points[point] == commonpt): #means that point belongs to a pentaagon
                    s_facet = self.ch.simplices[i]
                    facet_nlist, neigh_ind = self.facet_ncheck( facet_index =i, common_pt =commonpt)
                    while len(neigh_ind) < 7:
                        if len(neigh_ind) == 6:
                            facet_nlist.append(s_facet)
                            return facet_nlist, neigh_ind
                        neigh_facets, neigh_indices = self.facet_ncheck( facet_index =neigh_ind[-1], common_pt =commonpt )
                        neigh_ind = self.remove_duplicates(self.merge(list(neigh_ind),list(neigh_indices)))
                    return facet_nlist, neigh_ind
    
    def next_row_tri_v(self,  tri_index, tri_ind_cache):
        
        idx = range(len(self.ch.simplices))
        neighbours = self.ch.neighbors[tri_index]
        flat_neigh = neighbours.flatten()
        diff_list = []
        n_dict = []#{}
        for indx in tri_index:
            diff = list(set(self.ch.neighbors[indx]) - set(tri_ind_cache))
            diff_list.append(diff[0])
            for i in range(len(diff)):
                n_dict.append((indx,diff[i]))
        
        common = set(list(flat_neigh)) - set(list(tri_ind_cache))#- set(list(tri_index))

        return list(common), n_dict #vert_connections.reshape((-1))
        
        
    def next_row_tri_h(self,  tri_index, tri_ind_cache):
        
        idx = range(len(self.ch.simplices))
        neighbours = self.ch.neighbors[tri_index]
        flat_neigh = neighbours.flatten()
        diff_list = []
        x_vconn = []
        n_dict = []#{}
        for indx in tri_index:
            diff = list(set(self.ch.neighbors[indx]) - set(tri_ind_cache)) #list of faces left 
            
            diff_list.append(diff[0])
            for i in range(len(diff)):
                n_dict.append((indx,diff[i]))
                mid = list( (set(self.ch.neighbors[diff[i]]) & set(tri_index)) )#-set(indx))
                if len(mid) == 2:             
                    adj_diff = list(set(mid) - set([indx])) #list( (set(self.ch.neighbors[diff[i]]) & set(tri_index)) -set(indx))
                    pre_adj_diff = list(set(self.ch.neighbors[adj_diff[0]]) & set(tri_ind_cache))
                    print("preadjdiff", pre_adj_diff, adj_diff, mid, indx)
                    pre_ind = list(set(self.ch.neighbors[indx]) & set(tri_ind_cache))
                    print("preind", pre_ind,self.ch.neighbors[indx], tri_index)#_cache )
                    pre_vert = list(set(self.ch.neighbors[pre_adj_diff[0]]) & set(self.ch.neighbors[pre_ind[0]]))
                    print("prevert", pre_vert)
                    if len(pre_vert) != 0: #must have hit a pentagon
                        x_vconn.append((pre_vert[0], diff[i] ))
                else:
                    None

        common = set(list(flat_neigh)) - set(list(tri_ind_cache))#- set(list(tri_index))
        return list(common),n_dict, x_vconn #h_connections.reshape((-1))

        

    def construct_icosphere(self):
        self.icosahedron_vertices = self.normalize_v3(np.array(self.icosahedron_vertices))
        self.icosahedron_vertices.tolist()
        face_cache = []
        for face in self.faces:
            new_face_vert = [self.icosahedron_vertices[face[0]], self.icosahedron_vertices[face[1]], self.icosahedron_vertices[face[2]]]
            face_cache.append(new_face_vert)
        self.faces = np.asarray(face_cache )
    
        sph_vert, sph_tri = self.create_unit_sphere_vert()
        
        x, y, z = sph_vert[:,0], sph_vert[:,1], sph_vert[:,2] #Arrays corresponting the the cartesian coords
        
        self.ch = ConvexHull(sph_vert)
        self.ch = ConvexHull(self.ch.points[self.ch.vertices])
        return self.ch

    def plot_sphere(self, colours):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        x, y, z   = self.ch.points[self.ch.vertices][:,0],self.ch.points[self.ch.vertices][:,1], self.ch.points[self.ch.vertices][:,2]
        surf = ax.plot_trisurf(x,y,z, triangles=self.ch.simplices, cmap=plt.cm.Greys_r)
        surf.set_array(colours)
        self.icosahedron_vertices=np.asarray(self.icosahedron_vertices)
        ax.scatter(self.icosahedron_vertices[:,0],self.icosahedron_vertices[:,1], self.icosahedron_vertices[:,2], c='red')
        
        plt.show()
        return surf
    """
            
        if self.plot or self.replot:
            fig1 = plt.figure()
            fig2 = plt.figure() 
            
            
            self.ax = fig2.add_subplot(111, projection='3d')
            self.ax.set_title('3DLattice, Nu = %s' %(self.heart.p_fibrosis))
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')
            self.ax.set_zlabel('z')

           
            facecol = colours
            self.surf = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, facecolors=cm.gray(facecol), antialiased=False)
            #plt.draw()
            
            self.interval=25
            self.counter=0
            if animate == True:
                self.anim1 = animation.FuncAnimation(fig2, self.updatefig,
                            frames=200, interval=self.interval, blit=False)



    def updatefig(self, *args): #Function that yields the data for animation
   
        if (self.heart.time % self.heart.heartbeatsteps)==0 and self.heart.time!=0:    #why self.heart.time != 0??
            self.heart.excitecolumn()

        if self.plot==True and self.store==False:
            self.heart.onestep()
            
        if self.plot==True and self.store==True:
            self.gridintime.append(self.heart.grid)
            self.heart.onestep()
            
        self.ax.clear()
        self.im.set_array(self.heart.grid)
        facecol = np.absolute(self.heart.grid/self.heart.excitation)
        self.surf = self.ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,linewidth = 0, facecolors=cm.gray(facecol), antialiased=False)

        return self.surf,
        """
        


"""
s = Sphere(vertices = icosahedron_vertices, faces = faces, recursion_level = 4 )

s.construct_icosphere()
x, y, z   = s.ch.points[s.ch.vertices][:,0],s.ch.points[s.ch.vertices][:,1], s.ch.points[s.ch.vertices][:,2]


vertex1, vertex2, vertex3 = s.ch.points[s.ch.simplices[:,0]],s.ch.points[s.ch.simplices[:,1]], s.ch.points[s.ch.simplices[:,2]]# hull_pts[ch.simplices[:,1]], hull_pts[ch.simplices[:,2]]



#CHANGE TO SPHERICAL POLARS
sph_v_1 = s.sph_polar_convert(vertex1)
sph_v_2 = s.sph_polar_convert(vertex2) 
sph_v_3 = s.sph_polar_convert(vertex3) 

phi_avg = (sph_v_1[:,2] + sph_v_2[:,2] + sph_v_3[:,2])/3.
"""
"""            
pent_faces, pent_ind = s.find_pentagon()

val = 1
phi_avg[pent_ind] = val
phi_avg[phi_avg !=2] = 0

next_tri, vconn = s.next_row_tri_v(pent_ind, pent_ind) #Defines vertical connections 
phi_avg[next_tri] = val + 1
face_cache = np.hstack((pent_ind, next_tri)) 
next_tri2, hconn,  =  s.next_row_tri_h(next_tri, face_cache) 
phi_avg[next_tri2] = val + 1
for i in range(21):
    
    face_cache = np.hstack((face_cache, next_tri2))
    next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
    phi_avg[next_tri3] = val + i + 2
    face_cache = np.hstack((face_cache, next_tri3))
    next_tri4, hconn2 =  s.next_row_tri_h(next_tri3, face_cache)
    phi_avg[next_tri4] = val + i + 2
    next_tri2 = next_tri4
    vconn = vconn + vconn2#np.hstack((vconn, vconn2))
    #vconn = {**vconn, **vconn2}
    #hconn = {**hconn, **hconn2}
    hconn = hconn + hconn2#np.hstack((hconn, hconn2))

  
"""
"""
colours = phi_avg#sph_v_3[:,2] #Making the colors of the faces correspond to the value of the angle
colours = phi_avg/sum(phi_avg)
#rang =range(len(phi_avg));colours = np.array(rang)/sum(rang)
l =[513, 388, 773, 390, 1159, 774, 265, 778, 383, 780, 142, 146, 663, 24, 793, 1051, 290, 41, 45, 561, 562, 565, 315, 1086, 318, 321, 834, 835, 582, 470, 217, 97, 484, 486, 747, 110, 112, 507, 254, 763]
#print("vconn", vconn)
#print("hconn", hconn)
colours[l] = 1
colours[colours != 1] = 0
colours = colours/sum(colours)
suf =s.plot_sphere(colours)   

"""
    
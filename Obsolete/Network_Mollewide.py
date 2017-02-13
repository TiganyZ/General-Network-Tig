import numpy as np
from scipy.sparse import csr_matrix
import operator
from itertools import chain
import time
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import types
import pickle


def decision(p):
    if np.random.random()<p:
        return True
    else:
        return False
def trial():
    nodes=range(40000)
    a=np.arange(40000)
    b=np.arange(40000)
    c=np.arange(40000)
    np.random.shuffle(b)
    np.random.shuffle(c)
    l=[]
    for e in np.arange(40000):
        l.append((a[e],b[e]))
        l.append((c[e],b[e]))
        l.append((a[e],c[e]))
    array_transv=[]
    p=[0,1]
    return nodes,l,array_transv,1,p,0.5,0.5
        
    

class create_network:
    def __init__(self,array_nodesindices,array_vertical,array_transv,p_transv,impulse_start,p_dysf,p_unexcitable, excitation, hbs):
        """
        this inputs array_nodesindices as a list of faces indices,
        
        array_vertical and array_transv as  list of tuples
        """
        t0=time.time()
            
        # watch it because here we assume that faces go from 0 to n
        
        self.p_dysf=p_dysf
        self.p_transv=p_transv
        self.p_unexcitable=p_unexcitable
        self.excitation=excitation
        self.heartbeatssteps=hbs
        self.totalruns = 0
    
        self.array_vertical=array_vertical
        self.impulse_start=impulse_start
        self.size=len(array_nodesindices)
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]
        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
        self.array_alltransv=array_transv
               
        for elements in self.array_alltransv: #append transversal connections to a new list according to the probability of transv connect
          if decision(self.p_transv):
              self.array_transv.append(elements)
        
        self.connections = {}   #create dictionary with connections ex d={1:[0],2:[3,4],.....}
        for key,value in (self.array_vertical+self.array_transv):
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        
        for value,key in (self.array_vertical+self.array_transv): #append also connections from second to first elements
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        t1=time.time()
        print ('initialisation time:',(t1-t0))
     
    def excite(self,nodeindex):
        

     
        self.excited.append(nodeindex)
        self.nodes[nodeindex]=self.excitation

    def onestep(self):

        

        newnodes=np.zeros(self.size)
        time1=time.time()
    
        
        
        if len(self.excited)!=0:
            f = operator.itemgetter(*self.excited) # works well with list or array
            time2=time.time()
            #print("time taken",(time2-time1))
            
            
            self.excited=f(self.connections) # alternative function: list(self.connections[_] for _ in self.excited) check which one is faster
            try:
                self.excited=list(chain(*self.excited))#the output of f() is ( [....],[...],....) this changes into a unique list
            except TypeError:
                self.excited=list(chain(self.excited))
            time3=time.time()
            #print("time taken",(time3-time2)) #seems quick but can improve
        
            newnodes[self.excited]+=1   #this seems to take more than twice as much the previous operation 0.002
            #need to account for repeated elements
            time4=time.time()
            #print("time taken",(time4-time3))
        
            self.unexcited=np.random.rand(self.size)  #create array of dysf cells maybe can find a quicker way #this seems to take long 0.001
            self.unexcited[self.unexcited < self.p_unexcitable] = 1
            self.unexcited[self.unexcited != 1] = 0
            time5=time.time()
            #("time taken",(time5-time4))
            
            newnodes[newnodes>0]=1
            newnodes=newnodes-(self.dysf*self.unexcited*newnodes) #remove excited dysf cells which  are not excited
            newnodes*= (self.nodes==0) #removes refractory cells
        
            self.excited=np.flatnonzero(newnodes)
            self.excited=self.excited.tolist()
    
        time6=time.time()
        self.nodes-=1  
        self.nodes[self.nodes==-1]=0
        self.nodes += newnodes*self.excitation
        #print (np.flatnonzero(self.nodes))
        time7=time.time()
        #print("time taken",(time7-time6))
        
        #print ("total time", ( time7-time1))
   
    def reinitialise(self):
        self.nodes=np.zeros(self.size)
       
        self.array_transv=[]
        self.excited=[]
        self.dysf=np.random.rand(self.size)  #create array of dysf cells
        self.dysf[self.dysf < self.p_dysf] = 1
        self.dysf[self.dysf != 1] = 0
        self.unexcited=np.random.rand(self.size)
                
        for elements in self.array_alltransv: #append transversal connections to a new list according to the probability of transv connect
          if decision(self.p_transv):
              self.array_transv.append(elements)
         
        self.connections = {}   #create dictionary with connections ex d={1:[0],2:[3,4],.....}
        for key,value in (self.array_vertical+self.array_transv):
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
        
        for value,key in (self.array_vertical+self.array_transv): #append also connections from second to first elements
            try:
                self.connections[key].append(value)
            except KeyError:
                self.connections[key]=[value]
      
        

class run:
    
    def __init__(self,network, plot=False,store=True,runs=5000):
        
        
        self.runs=runs   
        self.plot= plot
        self.store=store
        self.network=network
        self.sph = []
        self.time=0
        self.nodeshistory=[]
        self.excitedhistory=[]
        self.nodeshistory.append(self.network.nodes)
        
        
    def propagate_n(self):
        for times in range(runs):
            if self.store==True:
                
                if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
                    for elements in self.network.impulse_start:
                        self.network.excite(elements)
                    
                self.time+=1
                self.network.totalruns+=1
                self.excitedhistory.append(self.network.excited)
            self.network.onestep()
            
    def propagate_a(self):
        for times in range(runs):
            if self.store==True:
                
                if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
                    for elements in self.network.impulse_start:
                        self.network.excite(elements)
                    
                self.time+=1
                self.network.totalruns+=1
                self.excitedhistory.append(self.network.excited)
            self.network.onestep()
            colours = self.network.nodes
            colours = colours/sum(colours) 
            yield colours
    
    def animator(self,sph):
        self.sph = sph
        if self.plot == True:
            self.surf, self.fig = self.plot_sphere_a(sph)
            self.sph = sph
            colours = self.propagate_a()
            self.anim1 = animation.FuncAnimation(self.fig, self.updatefig, fargs = (colours, self.surf),
                    frames=self.runs, interval=25, blit=False)
            plt.show()
          
            
                
                   
                
                
    def plot_sphere_a(self, sph):
        self.fig = plt.figure()
        #self.ax = self.fig.add_subplot(111, projection='3d')#.gca(projection='3d')projection = 'mollweide'
        self.ax = self.fig.add_subplot(111, projection = '3d') #'mollweide')#.gca(projection='3d')projection = 'mollweide'
        """self.ax.view_init(elev=20., azim=190)
        #self.ax.axis([-1,1,-1,1, -1, 1])
        
        self.ax.set_xlim(-0.55, 0.55)
        self.ax.set_ylim(-0.55, 0.55)
        self.ax.set_zlim(-0.55, 0.55)
        """
        ########
        #Maybe Project the 3d vertices to 2d
        ######
        #self.x, self.y, self.z   = sph.ch.points[sph.ch.vertices][:,0],sph.ch.points[sph.ch.vertices][:,1], sph.ch.points[sph.ch.vertices][:,2]
        self.x, self.y = sph.Mollewide_Projection()
        self.z = np.zeros(self.x.shape)
        
        self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=sph.ch.simplices, cmap=plt.cm.Greys_r, alpha = 1)
        self.triangles = sph.ch.simplices
        sph.icosahedron_vertices=np.asarray(sph.icosahedron_vertices)
        #self.ax.scatter(sph.icosahedron_vertices[:,0],sph.icosahedron_vertices[:,1], sph.icosahedron_vertices[:,2], c='red')
        
        plt.show()
        return self.surf, self.fig
        

    def updatefig(self, *args): #Function that yields the data for animation
        
        if self.network.totalruns%self.network.heartbeatssteps == 0: #self.time%self.network.heartbeatssteps==0:
            for elements in self.network.impulse_start:
                self.network.excite(elements)
                
        self.time+=1
        self.network.totalruns+=1
        self.excitedhistory.append(self.network.excited)
        self.network.onestep()
        colours = self.network.nodes
        colours = colours/sum(colours) 
        
        
        self.ax.clear()
        self.surf = self.ax.plot_trisurf(self.x,self.y,self.z, triangles=self.triangles, cmap=plt.cm.Greys_r, antialiased=False)
        self.surf.set_array(colours)
        """
        self.ax.set_xlim(-0.55, 0.55)
        self.ax.set_ylim(-0.55, 0.55)
        self.ax.set_zlim(-0.55, 0.55)
        """
        print("time", self.time)
        return self.surf,
        
                
import ico_propagation3 as sp


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

def num_func(recursion_level):
    t_0 = 3
    for x in range(recursion_level-2):
        t_n = 2*t_0 + 3
        t_0 = t_n
    return t_n
   
def define_connections(s):
    
    s.construct_icosphere()
    x, y, z   = s.ch.points[s.ch.vertices][:,0],s.ch.points[s.ch.vertices][:,1], s.ch.points[s.ch.vertices][:,2]
    vertex1 = s.ch.points[s.ch.simplices[:,0]]
    colours = vertex1
    pent_faces, pent_ind = s.find_pentagon(ind = 0)
    #pent_faces2, pent_ind2 = s.find_pentagon(ind = 3)
    next_tri, vconn = s.next_row_tri_v(pent_ind, pent_ind) #Defines vertical connections 
    face_cache = np.hstack((pent_ind, next_tri)) 
    next_tri2, hconn, x_v_conn=  s.next_row_tri_h(next_tri, face_cache) 
    vconn = vconn + x_v_conn
    """
    face_cache = np.hstack((face_cache, next_tri2))
    next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
    face_cache = np.hstack((face_cache, next_tri3))
    next_tri4, hconn2 , x_v_conn=  s.next_row_tri_h(next_tri3, face_cache)
    next_tri2 = next_tri4
    vconn = vconn + vconn2 + x_v_conn
    hconn = hconn + hconn2
    """
    """
    face_cache = np.hstack((face_cache, next_tri2))
    next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
    face_cache = np.hstack((face_cache, next_tri3))
    next_tri4, hconn2 , x_v_conn=  s.next_row_tri_h(next_tri3, face_cache)
    next_tri2 = next_tri4
    vconn = vconn + vconn2 + x_v_conn
    hconn = hconn + hconn2
    """
    startimp = np.hstack((face_cache, next_tri2))
       
    rang_val = int((num_func(s.recursion_level)+1)/2.)  #-2
    for i in range(rang_val):
        face_cache = np.hstack((face_cache, next_tri2))
        next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
        face_cache = np.hstack((face_cache, next_tri3))
        next_tri4, hconn2 , x_v_conn=  s.next_row_tri_h(next_tri3, face_cache)
        next_tri2 = next_tri4
        vconn = vconn + vconn2 + x_v_conn
        hconn = hconn + hconn2
    
    face_cache = np.hstack((face_cache, next_tri2))
    next_tri3, vconn2 =  s.next_row_tri_v(next_tri2, face_cache)
    vconn = vconn + vconn2
    
    
    
    pent_faces2, pent_ind2 = s.find_pentagon(ind = 3)
    next_trif, vconnf = s.next_row_tri_v(pent_ind2, pent_ind2) #Defines vertical connections 
    face_cachef = np.hstack((pent_ind2, next_trif)) 
    next_tri2f, hconnf, x_v_connf=  s.next_row_tri_h(next_trif, face_cachef) 
    vconnf = vconnf + x_v_connf
       
    rang_val = int((num_func(s.recursion_level)+1)/2.)
    for i in range(rang_val):
        face_cachef = np.hstack((face_cachef, next_tri2f))
        next_tri3f, vconn2f =  s.next_row_tri_v(next_tri2f, face_cachef)
        face_cachef = np.hstack((face_cachef, next_tri3f))
        next_tri4f, hconn2f , x_v_connf=  s.next_row_tri_h(next_tri3f, face_cachef)
        next_tri2f = next_tri4f
        vconnf = vconnf + vconn2f + x_v_connf
        hconnf = hconnf + hconn2f
    
    face_cachef = np.hstack((face_cachef, next_tri2f))
    next_tri3f, vconn2f =  s.next_row_tri_v(next_tri2f, face_cachef)
    vconnf = vconnf + vconn2f
    
    vconn = vconn + vconnf
    hconn = hconn + hconnf
    vconn = s.remove_duplicates(vconn)
    hconn = s.remove_duplicates(hconn)
    
    return colours, vconn, hconn, startimp#pent_ind


def connectome(s):
    vertex1, vertex2, vertex3 = s.ch.points[s.ch.simplices[:,0]],s.ch.points[s.ch.simplices[:,1]], s.ch.points[s.ch.simplices[:,2]]
    phi_avg = (vertex1+vertex2+vertex3)/3.
    
    linesv = [(tuple(phi_avg[i[0]]), tuple(phi_avg[i[1]]) )for i in vconn]
    colv = [(1,0,0,1) for i in range(len(linesv))]
    linesh = [(tuple(phi_avg[i[0]]), tuple(phi_avg[i[1]])) for i in hconn]
    colh = [(0,0,1,1) for i in range(len(linesh))]
    colx = [(0,1,0,1) for i in range(len(linesh))]
    coly = [(1,1,0,1) for i in range(len(linesh))]
    
    lines = linesv+linesh
    c = np.array(colv+colh)

    
    
    lc1 = Line3DCollection(linesv[int(len(linesv)/2):], colors=colv, linewidths=3)#[:int(5*len(linesv)/7)] [:int(len(linesv)/2)] + linesv[int(3*len(linesv)/4):]
    lc2 = Line3DCollection(linesh[int(len(linesh)/2):], colors=colh, linewidths=3)#[:int(5*len(linesh)/7)]
    #lc3 = Line3DCollection(linesv[int(len(linesv)/2):], colors=colx, linewidths=3)
    #lc4 = Line3DCollection(linesh[int(len(linesh)/2):], colors=coly, linewidths=3)
    #fig, ax = pl.subplots()
    #ax.add_collection(lc)
    fig2 = plt.figure()
    ax2 = fig2.gca(projection='3d')
    #ax2.plot(phi_avg[:,0], phi_avg[:,1], phi_avg[:,2], 'r', markersize=1)
    ax2.add_collection(lc1)
    ax2.add_collection(lc2)



s = sp.Sphere(vertices = icosahedron_vertices, faces = faces, recursion_level = 5 )
#
"""
s.construct_icosphere()

d=open('horiz_conn_rec_6xconn.pkl', 'rb')
hconn=pickle.load(d)

f=open('vert_conn_rec_6xconn.pkl', 'rb')
vconn=pickle.load(f)

e=open('startimp_ind_rec_6xconn.pkl', 'rb')
pent_ind=pickle.load(e)

g=open('colours_rec_6xconn.pkl', 'rb')
colours=pickle.load(g)

#h=open('sph_rec_6xconn.pkl', 'rb')
#s=pickle.load(h)

"""

colours, vconn, hconn, pent_ind = define_connections(s)

#fileobj1 = open('sph_rec_6.pkl', 'wb')
#pickle.dump(s, fileobj1, -1)
#fileobj1.close()
"""
fileobj1 = open('horiz_conn_rec_6xconn.pkl', 'wb')
pickle.dump(hconn, fileobj1, -1)

fileobj2 = open('vert_conn_rec_6xconn.pkl', 'wb')
pickle.dump(vconn, fileobj2, -1)

fileobj0 = open('startimp_ind_rec_6xconn.pkl', 'wb')
pickle.dump(pent_ind, fileobj0, -1)

fileobj3 = open('colours_rec_6xconn.pkl', 'wb')
pickle.dump(colours, fileobj3, -1)

fileobjsph = open('sph_rec_6xconn.pkl', 'wb')
pickle.dump(s, fileobjsph, -1)


fileobj1.close()
fileobj2.close()
fileobj0.close()
fileobj3.close()
fileobjsph.close()
"""


n = create_network(array_nodesindices = np.arange(len(colours)),
                   array_vertical = vconn,
                   array_transv = hconn,
                   p_transv = 0.15,
                   impulse_start = pent_ind,
                   p_dysf = 0.05,
                   p_unexcitable = 0.05,
                   excitation = 25, 
                   hbs = 110)
runc = run(network = n, plot=True,store=False,runs=1000)



#connectome(s)
#surface, figure = runc.plot_sphere_a(s)


#ax2.add_collection(lc3)
#ax2.add_collection(lc4)

#runc.plot(phi_avg[:,0], phi_avg[:,1], phi_avg[:,2], '.r', markersize=1)
#runc.ax.add_collection(lc1)
#runc.ax.add_collection(lc2)
#runc.ax.scatter(phi_avg[:,0], phi_avg[:,1], phi_avg[:,2], c='green')

runc.animator(s)

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=1800)
#runc.anim1.save('Recursion6xvom.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
#runc.anim1.save('')
#plt.show()

"""
d.close()
e.close()
f.close()
g.close()
#h.close()
"""
"""
phi = np.linspace(0, 2 * np.pi, 100)
theta = np.linspace(0, np.pi, 100)
xm =  0.9*np.outer(np.cos(phi), np.sin(theta)) 
ym =  0.9*np.outer(np.sin(phi), np.sin(theta))
zm =  0.9*np.outer(np.ones(np.size(phi)), np.cos(theta))
ax2.plot_surface(xm, ym, zm, cmap=plt.cm.Greys_r)
"""




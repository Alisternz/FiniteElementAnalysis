import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt

# Supporting Functions
def Psi1(x, L):
    return (1-x/L)
def Psi2(x, L):
    return (x/L)
def N1(x,L):
    return (1 - ((3*x**2) / L**2) + ((2*x**3) / L**3))
def N2(x,L):
    return ((x**3/L**2) - (2*x**2 / L) + x)
def N3(x,L):
    return ((3*x**2 / L**2) - (2*x**3 / L**3))
def N4(x,L):
    return ((x**3/L**2) - (x**2 / L))


# Info on System 
Shape1 = np.array([[0,0,0,1,0,0],
                   [0,0,0,0,1,0],
                   [0,0,0,0,0,1],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0,0]])

Shape2 = np.array([[1,0,0,0,0,0],
                   [0,1,0,0,0,0],
                   [0,0,1,0,0,0],
                   [0,0,0,1,0,0],
                   [0,0,0,0,1,0],
                   [0,0,0,0,0,1]])

Shape3 = np.array([[0,0,0,0,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0,0,0,0],
                   [1,0,0,0,0,0],
                   [0,1,0,0,0,0],
                   [0,0,1,0,0,0]])

Shapes = [Shape1,Shape2,Shape3]
Areas = [5e-4,5e-4,5e-4]
Es = [200e9,200e9,200e9]
Is = [1e-5,1e-5, 1e-5]

# Forces


# Uniform Distributed Load (UDL)
# Beam in question:


# Linearly Varying Load (LVL)



# Mid-Span Load




# Create Class
class Element:
    def __init__(self, E, A, I, Shape):
        self.E = E
        self.A = A
        self.I = I
        self.Shape = Shape

    def setNodes(self, x1, y1, L, angleDeg):
        self.x1 = x1
        self.y1 = y1
        self.angleDeg = angleDeg
        self.L = L
        self.angleRad = math.radians(self.angleDeg)
        self.x2 = self.x1 + self.L * np.cos(self.angleRad)
        self.y2 = self.y1 + self.L * np.sin(self.angleRad)


    def getTransformationMatrix(self):
        self.DOF = int(input("How many Degrees of Freedom; "))
        if (self.DOF == 4):
            k = (self.E * self.A) / self.L
            self.transformationMatrix = np.array([[1, -1], [-1, 1]])
        elif (self.DOF == 6):
            k = (self.E * self.I) / (self.L ** 3)
            B = (self.A * self.L**2) / self.I
            l = self.L
            self.transformationMatrix = np.array([[B,0,0,-B,0,0],
                                         [0, 12, 6*l, 0, -12, 6*l],
                                         [0, 6*l, 4 * l **2 , 0, -6*l, 2*l**2],
                                         [-B, 0, 0, B, 0, 0],
                                         [0, -12, -6*l, 0, 12, -6*l],
                                         [0, 6*l, 2*l**2, 0, -6*l, 4*l**2]])
        self.K = k * self.transformationMatrix

    def getkHat(self):
        angleRad = self.angleRad
        if (self.DOF == 4):
            self.lamda = np.array([[np.cos(angleRad), np.sin(angleRad), 0, 0],
                                         [0, 0, np.cos(angleRad), np.sin(angleRad)]])
        else:
            self.lamda = np.array([[np.cos(angleRad), np.sin(angleRad), 0, 0, 0, 0],
                                         [-np.sin(angleRad), np.cos(angleRad), 0, 0, 0, 0],
                                         [0, 0, 1, 0, 0, 0],
                                         [0,0,0, np.cos(angleRad), np.sin(angleRad), 0],
                                         [0, 0, 0, -np.sin(angleRad), np.cos(angleRad), 0],
                                         [0, 0, 0, 0, 0, 1]])
            
        self.kHat = self.lamda.T @ self.K @ self.lamda
        #print(self.kHat)

    def getkGE(self):
        self.kGe = self.Shape @ self.kHat @ self.Shape.T

    def updateValues(self, F, D, d):
        self.F = F
        self.D = D
        self.d = d
    
    def transverseDeflection(self, x, d):
        transverse = d[1][0]*N1(x, self.L) + d[2][0]*N2(x, self.L) + d[4][0]*N3(x, self.L) + d[5][0]*N4(x, self.L)
        return transverse
    
    def axialDeflection(self, x, d):
        axial = d[0][0]*Psi1(x, self.L) + d[3][0]*Psi2(x, self.L)
        return axial

    def totalDeflection(self, transverse, axial):
        xDeflection = axial * np.cos(self.angleRad) - transverse * np.sin(self.angleRad)
        yDeflection = axial * np.sin(self.angleRad) + transverse * np.cos(self.angleRad)
        return xDeflection,yDeflection

    def getQUDL(self, w):
        self.f_UDL = np.array([[0],
                      [(w*self.L)/2],
                      [(w*self.L**2)/12],
                      [0],
                      [(w*self.L)/2],
                      [(-w*self.L**2)/12]])
        self.F_UDL = self.lamda.T @ self.f_UDL
        self.Q_UDL = self.Shape @ self.F_UDL
    
    def getQLVL(self, w):
        self.f_LVL = np.array([[0],
                               [3*w*self.L / 20],
                               [w*self.L**2 / 30],
                               [0],
                               [7*w*self.L / 20],
                               [-w*self.L**2 / 20]])
        self.F_LVL = self.lamda.T @ self.f_LVL
        self.Q_LVL = self.Shape @ self.F_LVL

    def getQMSL(self, w):
        self.f_MSL = np.array([[0],
                               [w/2],
                               [w*self.L/8],
                               [0],
                               [w/2],
                               [-w*self.L/8]])
        self.F_MSL = self.lamda.T @ self.f_MSL
        self.Q_MSL = self.Shape @ self.F_MSL


        
        
        

# Initialise Elements
numElements = int(input("Number of Elements: "))
Elements = [Element(Es[i], Areas[i], Is[i], Shapes[i]) for i in range(numElements)]


# Element: (x,y), Length, Angle(deg)
Elements[0].setNodes(0, 0, 3, 90)
Elements[1].setNodes(0, 3, 4.5, 0)
Elements[2].setNodes(4.5, 3, 3, -90)


# Fill Elements (BASED ON NODES and LOADING)
for i in range(numElements):
    Elements[i].getTransformationMatrix()
    Elements[i].getkHat()
    Elements[i].getkGE()
    if (i==0):
        overallKG = Elements[i].kGe
        #print("Create KG")
    else:
        overallKG += Elements[i].kGe
        #print("Add to KG")
    
    
    
# Calculate Q from Nodal, UDL, MSL, LVL

# Loads

Elements[0].getQUDL(0)
Elements[0].getQLVL(0)
Elements[0].getQMSL(0)

Elements[1].getQUDL(-10e3)
Elements[1].getQLVL(0)
Elements[1].getQMSL(-50e3)

Elements[2].getQUDL(0)
Elements[2].getQLVL(0)
Elements[2].getQMSL(0)

Q = np.array([[10e3],
              [0],
              [0],
              [10e3],
              [0],
              [0]])

for i in range(numElements):
    Q = Q + Elements[i].Q_UDL + Elements[i].Q_LVL + Elements[i].Q_MSL


# Calculate overall Quantities
q = sp.linalg.solve(overallKG, Q)

# Get Info on Each Element
for i in range(numElements):
    F = Elements[i].kHat @ Elements[i].Shape.T @ q
    print(f"Element {i}: {F}")
    # print(f"Element {i} Angle(Degs): {Elements[i].angleDeg}")
    # print(f"Element {i} Angle(Rads): {Elements[i].angleRad}")
    # print(f"Element {i} K: {Elements[i].K}")
    # print(f"Element {i} kHat: {Elements[i].kHat}")
    # print(f"Element {i} KgE: {Elements[i].kGe}")
    #print(f"KG: {overallKG}")
    # print(f"Element {i} F: {F}")
    D = Elements[i].Shape.T @ q
    d = Elements[i].lamda @ D
    
    Elements[i].updateValues(F, D, d)
    print(f"{Elements[i].d} Deflection OF ELEMENT {i}")
# print(f"q: {q}")
#print(overallKG, Q)



# Plotting

# numPoints = 50
# for i in range(numElements):
#     coords =[[], #xs
#             []]  #ys
#     ls = np.linspace(0, Elements[i].L, numPoints) #Points along length of bar
#     for l in ls:
#         transverse = Elements[i].transverseDeflection(l,d)
#         axial = Elements[i].axialDeflection(l,d)
#         xDeflection, yDeflection = Elements[i].totalDeflection(transverse, axial)
#         if (l==3):
#           print(f"x= {coords[0]}, y={coords[1]}")
#           print(xDeflection)
#         coords[0].append(Elements[i].x1 + l*np.cos(Elements[i].angleRad) + 50 * xDeflection)
#         coords[1].append(Elements[i].y1 + l*np.sin(Elements[i].angleRad) + 50 * yDeflection)
#         if (l==3):
#           print(coords[0], coords[1])
            
#     plt.plot(coords[0], coords[1])
#     plt.grid(True)
#     #print()

N_points = 40
disp_multi = 30
for j in range(numElements):
    transverse = Elements[j].transverseDeflection(Elements[j].L ,d)
    axial = Elements[j].axialDeflection(Elements[j].L,d)
    print(Elements[j].d[1][0], "hi")
    xDeflection, yDeflection = Elements[j].totalDeflection(transverse, axial)
    node1XG = Elements[j].x1
    node1YG = Elements[j].y1
    node2XG = Elements[j].x2
    node2YG = Elements[j].y2
    x_e = np.linspace(0, Elements[j].L, N_points)

    ax_1 = (1-(x_e/Elements[j].L))

    ax_2 = (x_e/Elements[j].L)

    N_1 = 1 - (3*x_e**2)/(Elements[j].L**2) + (2*x_e**3)/(Elements[j].L**3)
    N_2 = (x_e**3)/(Elements[j].L**2) - (2*x_e**2)/(Elements[j].L) + x_e
    N_3 = (3*x_e**2)/(Elements[j].L**2) - (2*x_e**3)/(Elements[j].L**3)
    N_4 = (x_e**3)/(Elements[j].L**2) - (x_e**2)/Elements[j].L

    u_x = ax_1 * Elements[j].d[0][0] + ax_2 * Elements[j].d[3][0]
    v_x = N_1 * Elements[j].d[1][0] + N_2 * Elements[j].d[2][0] + N_3 * Elements[j].d[4][0] + N_4 * Elements[j].d[5][0]

    defl_XG = u_x * np.cos(np.radians(Elements[j].angleDeg)) - v_x * np.sin(np.radians(Elements[j].angleDeg))
    defl_YG = u_x * np.sin(np.radians(Elements[j].angleDeg)) + v_x * np.cos(np.radians(Elements[j].angleDeg))

    undefl_base_XG = np.linspace(node1XG, node2XG, N_points)
    undefl_base_YG = np.linspace(node1YG, node2YG, N_points)

    deflected_XG = undefl_base_XG + disp_multi * defl_XG
    deflected_YG = undefl_base_YG + disp_multi * defl_YG

    plt.plot(undefl_base_XG, undefl_base_YG, 'b--')
    plt.plot(deflected_XG, deflected_YG, 'r')
    
plt.show()






import numpy as np

def get_evals(M):
    """Function returning the sorted eigenvalues of the matrix M"""
    evals, evecs = np.linalg.eig(M)
    return sorted(evals)
def find_degenerates(energies):
    """Function returning a list of the energies of the Hückel orbitals and a list with the corresponding degeneracies"""
    degen_energies=[]
    degeneracies=[]
    degen=1
    for i in range(0,len(energies)-1):
        if not np.isclose(energies[i],energies[i+1],atol=0.00001):
            degen_energies.append(energies[i])
            degeneracies.append(degen)
            degen=1
        else:
            degen+=1
    degen_energies.append(energies[len(energies)-1])
    degeneracies.append(degen)
    return degen_energies, degeneracies
    
# Specification of problem:
print("Huckel MO calculator")
try:
    category=int(input("What system do you want to solve? Press 1: linear polyene, 2: cyclic polyene, 3: platonic solids, 4: Buckminsterfullerene "))
except:
    print("Please enter a valid category number.")
assert category>0 and category<5, "Please enter a valid category number."

#Linear polyene with n carbons
if category==1:
    try:
        n=int(input("Number of carbon atoms in the linear polyene: "))
    except:
        print("Bug in user input.")

    H=np.zeros((n,n))
    for i in range(0,n-1):
        H[i][i+1]=-1 # -1 if two atoms are adjacent 
    H=H+np.transpose(H)

    energies=get_evals(H)
    degen_energies, degeneracies = find_degenerates(energies)

    print("The energies of the orbitals are: "+str(degen_energies)+" with degeneracies: "+str(degeneracies))

#Cyclic polyene with n carbon
if category==2:
    try:
        n=int(input("Number of carbon atoms in the cyclic polyene: "))
    except:
        print("Bug in user input.")

    H=np.zeros((n,n))
    for i in range(0,n):
        H[i][(i+1)%n]=-1 #-1 if two atoms are adjacent modulo n
    H=H+np.transpose(H)
    #print(H)
    energies=get_evals(H)
    degen_energies, degeneracies = find_degenerates(energies)


    print("The energies of the orbitals are: "+str(degen_energies)+" with degeneracies: "+str(degeneracies))

#Platonic solids
if category==3:
    structure=str(input("Which Platonic solid are you interested in? (tetrahedron, octahedron, cube, icosahedron, dodecahedron) "))
    assert structure=="tetrahedron" or structure=="cube" or structure=="dodecahedron" or structure=="octahedron" or structure=="icosahedron", "Not a valid Platonic solid"

    if structure=="tetrahedron":
        H=np.zeros((4,4))
        for i in range(0,3):
            for j in range(i+1,4):
                H[i][j]=-1
        H=H+np.transpose(H)

    elif structure=="cube":
        H=np.zeros((8,8))
        for i in range(0,4):
            H[i][(i+1)%4]=-1      #first square
            H[i+4][(i+1)%4+4]=-1  #second square
            H[i][i+4]=-1          #connections between squares
        H=H+np.transpose(H)

    elif structure=="dodecahedron":
        H=np.zeros((20,20))
        for i in range(0,5):
            H[i][(i+1)%5]=-1         #inner circle
            H[i+15][(i+16)%5+15]=-1  #outer circle
            H[i][2*i+5]=-1           #inner and middle circle
            H[2*i+6][i+15]=-1        #middle and outer circle
        for i in range(0,10):
            H[i+5][(i+1)%10+5]=-1    #middle circle
        H=H+np.transpose(H)
    
    elif structure=="octahedron":
        H=np.ones((6,6))
        H=-1*H
        for i in range(0,6):
            H[i][i]=0
        H[0][4]=0
        H[1][5]=0
        H[2][3]=0
        H[4][0]=0
        H[5][1]=0
        H[3][2]=0
    
    elif structure=="icosahedron":
        H=np.zeros((12,12))
        for i in range(0,3):
            H[i][(i+1)%3]=-1            # first circle (in the middle)
            H[i+9][(i+1)%3+9]=-1        # third (outermost) circle 
            H[i][2*i+3]=-1              # first and second circle
            H[i][2*i+4]=-1              # first and second circle
            H[0][8]=-1                  # first and second circle
            H[1][4]=-1                  # first and second circle
            H[2][6]=-1                  # first and second circle
            H[i+9][2*i+3]=-1            # second and third circle
            H[i+9][2*i+4]=-1            # second and third circle
            H[9][5]=-1                  # second and third circle
            H[10][7]=-1                 # second and third circle
            H[11][3]=-1                 # second and third circle
        for i in range(0,6):
            H[i+3][(i+1)%6+3]=-1        # second circle 
        H=H+np.transpose(H)

    energies=get_evals(H)
    degen_energies, degeneracies = find_degenerates(energies)

    print("The energies of the orbitals are: "+str(degen_energies)+" with degeneracies: "+str(degeneracies))

#Buckminsterfullerene
if category==4:
    H=np.zeros((60,60))   
    for i in range(0,5):
        H[i][(i+1)%5]=-1            # first circle (in the middle)
        H[i+55][(i+1)%5+55]=-1      # fifth (outermost) circle
        H[i][3*i+5]=-1              # first and second circle
        H[3*i+6][4*i+20]=-1         # second and third circle
        H[3*i+7][4*i+23]=-1         # second and third circle
        H[4*i+21][3*i+40]=-1        # third and fourth circle
        H[4*i+22][3*i+42]=-1        # third and fourth circle
        H[3*i+41][i+55]=-1          # fourth and fifth circle
    for i in range(0,15):
        H[i+5][(i+1)%15+5]=-1       # second circle
        H[i+40][(i+1)%15+40]=-1     # fourth circle
    for i in range(0,20):
        H[i+20][(i+1)%20+20]=-1     #third circle

    H=H+np.transpose(H)
    energies=get_evals(H)
    degen_energies, degeneracies = find_degenerates(energies)

    print("The energies of the orbitals are: "+str(degen_energies)+" with degeneracies: "+str(degeneracies))


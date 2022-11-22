import numpy as np
import matplotlib.pyplot as plt
import math
import random

bvals=[100.0];
bnum=len(bvals)

L=10
D=3
usefultrials=15000
wastetrials=15000

V11=np.zeros(bnum);
V13=np.zeros(bnum);

I2=np.zeros((2,2), dtype=complex)
sig1=np.zeros((2,2), dtype=complex)
sig2=np.zeros((2,2), dtype=complex)
sig3=np.zeros((2,2), dtype=complex)
sig1[0,0]=0.0
sig1[0,1]=1.0
sig1[1,0]=1.0
sig1[1,1]=0.0
sig2[0,0]=0.0
sig2[0,1]=-1.0j
sig2[1,0]=1.0j
sig2[1,1]=0.0
sig3[0,0]=1.0
sig3[0,1]=0.0
sig3[1,0]=0.0
sig3[1,1]=-1.0
I2[0,0]=1.0
I2[0,1]=0.0
I2[1,0]=0.0
I2[1,1]=1.0


U=np.zeros((2,2),dtype=complex)
Udag=np.zeros((2,2),dtype=complex)
UUb=np.zeros((2,2),dtype=complex)
Ubdag=np.zeros((2,2),dtype=complex)
k=0.1+0.1j
Ub=np.zeros((2,2),dtype=complex)
U2=np.zeros((2,2),dtype=complex)
U3=np.zeros((2,2),dtype=complex)
U4=np.zeros((2,2),dtype=complex)
U5=np.zeros((2,2),dtype=complex)
U6=np.zeros((2,2),dtype=complex)
U7=np.zeros((2,2),dtype=complex)
U3dag=np.zeros((2,2),dtype=complex)
U4dag=np.zeros((2,2),dtype=complex)
U6dag=np.zeros((2,2),dtype=complex)
U7dag=np.zeros((2,2),dtype=complex)
config=np.zeros((L,L,L,D,4))


def sweep( numoftrials):
	V11ret=0
	V13ret=0
	for trial in range(numoftrials):
		for x in range(L):
			for y in range(L):
				for t in range(L):
					for mu in range(D):
						U=config[x,y,t,mu,0]*I2+1.0j*config[x,y,t,mu,1]*sig1+1.0j*config[x,y,t,mu,2]*sig2+1.0j*config[x,y,t,mu,3]*sig3

						Ub=np.zeros((2,2),dtype=complex)
						for nu in range(D):						
							if(nu==mu):
								continue;

							r=[x,y,t];
							rplusmu=[x,y,t];
							rplusnu=[x,y,t];
							rminnu=[x,y,t];
							rplmuminnu=[x,y,t];
							rplusmu[mu]=(r[mu]+1)%L
							rplusnu[nu]=(r[nu]+1)%L
							rminnu[nu]=(r[nu]-1)%L
							rplmuminnu[mu]=(r[mu]+1)%L
							rplmuminnu[nu]=(r[nu]-1)%L

							U2=config[rplusmu[0],rplusmu[1],rplusmu[2],mu,0]*I2+1.0j*config[rplusmu[0],rplusmu[1],rplusmu[2],nu,1]*sig1+1.0j*config[rplusmu[0],rplusmu[1],rplusmu[2],nu,2]*sig2+1.0j*config[rplusmu[0],rplusmu[1],rplusmu[2],nu,3]*sig3
							U3=config[rplusnu[0],rplusnu[1],rplusnu[2],mu,0]*I2+1.0j*config[rplusnu[0],rplusnu[1],rplusnu[2],mu,1]*sig1+1.0j*config[rplusnu[0],rplusnu[1],rplusnu[2],mu,2]*sig2+1.0j*config[rplusnu[0],rplusnu[1],rplusnu[2],mu,3]*sig3
							U4=config[r[0],r[1],r[2],nu,0]*I2+1.0j*config[r[0],r[1],r[2],nu,1]*sig1+1.0j*config[r[0],r[1],r[2],nu,2]*sig2+1.0j*config[r[0],r[1],r[2],nu,3]*sig3
							U5=config[rminnu[0],rminnu[1],rminnu[2],nu,0]*I2+1.0j*config[rminnu[0],rminnu[1],rminnu[2],nu,1]*sig1+1.0j*config[rminnu[0],rminnu[1],rminnu[2],nu,2]*sig2+1.0j*config[rminnu[0],rminnu[1],rminnu[2],nu,3]*sig3
							U6=config[rminnu[0],rminnu[1],rminnu[2],mu,0]*I2+1.0j*config[rminnu[0],rminnu[1],rminnu[2],mu,1]*sig1+1.0j*config[rminnu[0],rminnu[1],rminnu[2],mu,2]*sig2+1.0j*config[rminnu[0],rminnu[1],rminnu[2],mu,3]*sig3
							U7=config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu,0]*I2+1.0j*config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu,1]*sig1+1.0j*config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu,2]*sig2+1.0j*config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu,3]*sig3

							U3dag[0][0]=np.conj(U3[0][0])
							U3dag[0][1]=np.conj(U3[1][0])
							U3dag[1][0]=np.conj(U3[0][1])
							U3dag[1][1]=np.conj(U3[1][1])
							U4dag[0][0]=np.conj(U4[0][0])
							U4dag[0][1]=np.conj(U4[1][0])
							U4dag[1][0]=np.conj(U4[0][1])
							U4dag[1][1]=np.conj(U4[1][1])
							U6dag[0][0]=np.conj(U6[0][0])
							U6dag[0][1]=np.conj(U6[1][0])
							U6dag[1][0]=np.conj(U6[0][1])
							U6dag[1][1]=np.conj(U6[1][1])
							U7dag[0][0]=np.conj(U7[0][0])
							U7dag[0][1]=np.conj(U7[1][0])
							U7dag[1][0]=np.conj(U7[0][1])
							U7dag[1][1]=np.conj(U7[1][1])

							Ub=Ub+np.matmul(U2,np.matmul(U3dag,U4dag))
							Ub=Ub+np.matmul(U7dag,np.matmul(U6dag,U5))

						k=abs(math.sqrt(np.linalg.det(Ub)))
						Ub=Ub/k

						a0=0.00
						accepted=0;
						while(accepted==1):
							xval=random.uniform(math.exp(-4*b*k),1)
							a0=1.0+1.0/(2*b*k)*math.log(xval)

							if(random.random()<math.sqrt(1.0-a0*a0)):
								accepted=1

						a1=random.gauss(0, 1)
						a2=random.gauss(0, 1)
						a3=random.gauss(0, 1)

						maginv=math.sqrt(1-a0*a0)/math.sqrt(a1*a1+a2*a2+a3*a3)

						a1=maginv*a1
						a2=maginv*a2
						a3=maginv*a3

						UUb=a0*I2+1.0j*a1*sig1+1.0j*a2*sig2+1.0j*a3*sig3

						Ubdag[0][0]=np.conj(Ub[0][0])
						Ubdag[0][1]=np.conj(Ub[1][0])
						Ubdag[1][0]=np.conj(Ub[0][1])
						Ubdag[1][1]=np.conj(Ub[1][1])

						U=np.matmul(UUb,Ubdag)

						config[x,y,t,mu,0]=0.5*np.trace(U).real
						config[x,y,t,mu,1]=0.5*np.trace(np.matmul(U, -1.0j*sig1)).real
						config[x,y,t,mu,2]=0.5*np.trace(np.matmul(U, -1.0j*sig2)).real
						config[x,y,t,mu,3]=0.5*np.trace(np.matmul(U, -1.0j*sig3)).real

		U=config[5,5,5,1,0]*I2+1.0j*config[5,5,5,1,1]*sig1+1.0j*config[5,5,5,1,2]*sig2+1.0j*config[5,5,5,1,3]*sig3
		Udag[0][0]=np.conj(U[0][0])
		Udag[0][1]=np.conj(U[1][0])
		Udag[1][0]=np.conj(U[0][1])
		Udag[1][1]=np.conj(U[1][1])


		V11ret=V11ret+np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig1,Udag))))).real*np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig1,Udag))))).real;
		V13ret=V13ret+np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig3,Udag))))).real*np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig3,Udag))))).real;

	V11ret=V11ret/numoftrials
	V13ret=V13ret/numoftrials

	return V11ret, V13ret;

for runnum in range(bnum):
	b=bvals[runnum];															#b=-2/g2
	random.seed()
	
	for x in range(L):
		for y in range(L):
			for t in range(L):
				for mu in range(D):
					config[x,y,t,mu,0]=random.gauss(0.0,1.1)
					config[x,y,t,mu,1]=random.gauss(0.0,1.1)
					config[x,y,t,mu,2]=random.gauss(0.0,1.1)
					config[x,y,t,mu,3]=random.gauss(0.0,1.1)

					mag=np.sqrt(config[x,y,t,mu,0]*config[x,y,t,mu,0]+config[x,y,t,mu,1]*config[x,y,t,mu,1]+config[x,y,t,mu,2]*config[x,y,t,mu,2]+config[x,y,t,mu,3]*config[x,y,t,mu,3])

					config[x,y,t,mu,0]=config[x,y,t,mu,0]/mag
					config[x,y,t,mu,1]=config[x,y,t,mu,1]/mag
					config[x,y,t,mu,2]=config[x,y,t,mu,2]/mag
					config[x,y,t,mu,3]=config[x,y,t,mu,3]/mag					

	for x in range(L):
		for y in range(L):
			for t in range(L):
				for mu in range(D):
					config[x,y,t,mu,:]=1.0*config[x,y,t,mu,:]/math.sqrt(config[x,y,t,mu,0]*config[x,y,t,mu,0]+config[x,y,t,mu,1]*config[x,y,t,mu,1]+config[x,y,t,mu,2]*config[x,y,t,mu,2]+config[x,y,t,mu,3]*config[x,y,t,mu,3])

	sweep( wastetrials)

	V11[runnum], V13[runnum]=sweep( usefultrials)
	print(b,V11[runnum],V13[runnum])
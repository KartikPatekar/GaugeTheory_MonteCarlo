import numpy as np
import matplotlib.pyplot as plt
import math
import random
from progress.bar import Bar

#Convention of matrix: 0=theta , 1=psi , 3=phi

bvals=[1000,100, 10,1,0.05];
bnum=len(bvals)

dangle=0.1
L=10
D=3
usefultrials=8000
wastetrials=10000
bar = Bar('Processing', max=bnum*(usefultrials+wastetrials))

V11=np.zeros(bnum);
V13=np.zeros(bnum);

sig1=np.zeros((2,2), dtype=complex)
sig3=np.zeros((2,2), dtype=complex)
sig1=[[0,1],[1,0]]
sig3=[[1,0],[0,-1]]

for runnum in range(bnum):
	b=bvals[runnum];															#b=-2/g2
	random.seed()
	userandnums=np.random.rand(usefultrials);
	wasterandnums=np.random.rand(wastetrials);

	randth=math.pi*np.random.rand(L,L,L,D);
	randpsi=math.pi*np.random.rand(L,L,L,D);
	randphi=2*math.pi*np.random.rand(L,L,L,D);
	config=np.zeros((L,L,L,D,2,2), dtype=complex);
	config[:,:,:,:,0,0]=np.cos(randth)+1.0j*np.sin(randth)*np.sin(randpsi)*np.sin(randphi)
	config[:,:,:,:,0,1]=np.sin(randth)*np.cos(randpsi)+1.0j*np.sin(randth)*np.sin(randpsi)*np.cos(randphi)
	config[:,:,:,:,1,0]=-np.sin(randth)*np.cos(randpsi)+1.0j*np.sin(randth)*np.sin(randpsi)*np.cos(randphi)
	config[:,:,:,:,1,1]=np.cos(randth)-1.0j*np.sin(randth)*np.sin(randpsi)*np.sin(randphi)

	randth=math.pi*np.random.rand(500)*dangle;
	randpsi=math.pi*np.random.rand(500)*dangle;
	randphi=2*math.pi*np.random.rand(500)*dangle;
	Urandmatrices=np.zeros((1000,2,2), dtype=complex);
	Urandmatrices[0:500,0,0]=np.cos(randth)+1.0j*np.sin(randth)*np.sin(randpsi)*np.sin(randphi)
	Urandmatrices[0:500,0,1]=np.sin(randth)*np.cos(randpsi)+1.0j*np.sin(randth)*np.sin(randpsi)*np.cos(randphi)
	Urandmatrices[0:500,1,0]=-np.sin(randth)*np.cos(randpsi)+1.0j*np.sin(randth)*np.sin(randpsi)*np.cos(randphi)
	Urandmatrices[0:500,1,1]=np.cos(randth)-1.0j*np.sin(randth)*np.sin(randpsi)*np.sin(randphi)
	Urandmatrices[500:1000,0,0]=np.conj(Urandmatrices[0:500,0,0])
	Urandmatrices[500:1000,0,1]=np.conj(Urandmatrices[0:500,1,0])
	Urandmatrices[500:1000,1,0]=np.conj(Urandmatrices[0:500,0,1])
	Urandmatrices[500:1000,1,1]=np.conj(Urandmatrices[0:500,1,1])

	U1old=np.zeros((2,2), dtype=complex)
	U1new=np.zeros((2,2), dtype=complex)
	Urand=np.zeros((2,2), dtype=complex)
	U2=np.zeros((2,2), dtype=complex)
	U3=np.zeros((2,2), dtype=complex)
	U4=np.zeros((2,2), dtype=complex)
	U5=np.zeros((2,2), dtype=complex)
	U6=np.zeros((2,2), dtype=complex)
	U7=np.zeros((2,2), dtype=complex)
	U3dag=np.zeros((2,2), dtype=complex)
	U4dag=np.zeros((2,2), dtype=complex)
	U6dag=np.zeros((2,2), dtype=complex)
	U7dag=np.zeros((2,2), dtype=complex)
	U=np.zeros((2,2), dtype=complex)
	Udag=np.zeros((2,2), dtype=complex)

	randth=0;
	randpsi=0;
	randphi=0;
	r=[0,0,0];
	rplusmu=[0,0,0];
	rplusnu=[0,0,0];
	rminnu=[0,0,0];
	rplmuminnu=[0,0,0];
	Sold=0;
	Snew=0;

	for trial in range(wastetrials):
		bar.next()
		for x in range(L):
			for y in range(L):
				for t in range(L):
					for mu in range(D):
						# randth=math.pi*random.random()*dangle;
						# randpsi=math.pi*random.random()*dangle;
						# randphi=2*math.pi*random.random()*dangle;
						# Urand[0,0]=math.cos(randth)+1.0j*math.sin(randth)*math.sin(randpsi)*math.sin(randphi)
						# Urand[0,1]=math.sin(randth)*math.cos(randpsi)+1.0j*math.sin(randth)*math.sin(randpsi)*math.cos(randphi)
						# Urand[1,0]=-math.sin(randth)*math.cos(randpsi)+1.0j*math.sin(randth)*math.sin(randpsi)*math.cos(randphi)
						# Urand[1,1]=math.cos(randth)-1.0j*math.sin(randth)*math.sin(randpsi)*math.sin(randphi)
						Urand=Urandmatrices[random.randrange(1000)];

						U1old=config[x,y,t,mu];
						U1new=np.matmul(Urand, U1old);
						Sold=0
						Snew=0

						for nu in range(D):
							if(nu==mu):
								continue;

							r=[x,y,t];
							rplusmu=[x,y,t];
							rplusnu=[x,y,t];
							rminnu=[x,y,t];
							rplmuminnu==[x,y,t];
							rplusmu[mu]=(r[mu]+1)%L
							rplusnu[nu]=(r[nu]+1)%L
							rminnu[nu]=(r[nu]-1)%L
							rplmuminnu[mu]=(r[mu]+1)%L
							rplmuminnu[nu]=(r[nu]-1)%L

							U2=config[rplusmu[0],rplusmu[1],rplusmu[2],nu]
							U3=config[rplusnu[0],rplusnu[1],rplusnu[2],mu]
							U4=config[r[0],r[1],r[2],nu]
							U5=config[rminnu[0],rminnu[1],rminnu[2],nu]
							U6=config[rminnu[0],rminnu[1],rminnu[2],mu]
							U7=config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu]

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

							Sold=Sold-1.0*b*np.trace(np.matmul(U1old,np.matmul(U2,np.matmul(U3dag,U4dag)))).real
							Snew=Snew-1.0*b*np.trace(np.matmul(U1new,np.matmul(U2,np.matmul(U3dag,U4dag)))).real
							Sold=Sold-1.0*b*np.trace(np.matmul(U1old,np.matmul(U7dag,np.matmul(U6dag,U5)))).real
							Snew=Snew-1.0*b*np.trace(np.matmul(U1new,np.matmul(U7dag,np.matmul(U6dag,U5)))).real
						
						if(Sold<Snew):
							if(wasterandnums[trial]<=math.exp(Sold-Snew)):
								config[x,y,t,mu]=U1new;

						if(Sold>Snew):
								config[x,y,t,mu]=U1new;				

	for trial in range(usefultrials):
		bar.next()
		for x in range(L):
			for y in range(L):
				for t in range(L):
					for mu in range(D):
						randth=math.pi*random.random()*dangle;
						randpsi=math.pi*random.random()*dangle;
						randphi=2*math.pi*random.random()*dangle;
						Urand[0,0]=math.cos(randth)+1.0j*math.sin(randth)*math.sin(randpsi)*math.sin(randphi)
						Urand[0,1]=math.sin(randth)*math.cos(randpsi)+1.0j*math.sin(randth)*math.sin(randpsi)*math.cos(randphi)
						Urand[1,0]=-math.sin(randth)*math.cos(randpsi)+1.0j*math.sin(randth)*math.sin(randpsi)*math.cos(randphi)
						Urand[1,1]=math.cos(randth)-1.0j*math.sin(randth)*math.sin(randpsi)*math.sin(randphi)

						U1old=config[x,y,t,mu];
						U1new=np.matmul(Urand, U1old);
						Sold=0
						Snew=0

						for nu in range(D):
							if(nu==mu):
								continue;

							r=[x,y,t];
							rplusmu=[x,y,t];
							rplusnu=[x,y,t];
							rminnu=[x,y,t];
							rplmuminnu==[x,y,t];
							rplusmu[mu]=(r[mu]+1)%L
							rplusnu[nu]=(r[nu]+1)%L
							rminnu[nu]=(r[nu]-1)%L
							rplmuminnu[mu]=(r[mu]+1)%L
							rplmuminnu[nu]=(r[nu]-1)%L

							U2=config[rplusmu[0],rplusmu[1],rplusmu[2],nu]
							U3=config[rplusnu[0],rplusnu[1],rplusnu[2],mu]
							U4=config[r[0],r[1],r[2],nu]
							U5=config[rminnu[0],rminnu[1],rminnu[2],nu]
							U6=config[rminnu[0],rminnu[1],rminnu[2],mu]
							U7=config[rplmuminnu[0],rplmuminnu[1],rplmuminnu[2],nu]

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

							Sold=Sold-1.0*b*np.trace(np.matmul(U1old,np.matmul(U2,np.matmul(U3dag,U4dag)))).real
							Snew=Snew-1.0*b*np.trace(np.matmul(U1new,np.matmul(U2,np.matmul(U3dag,U4dag)))).real
							Sold=Sold-1.0*b*np.trace(np.matmul(U1old,np.matmul(U7dag,np.matmul(U6dag,U5)))).real
							Snew=Snew-1.0*b*np.trace(np.matmul(U1new,np.matmul(U7dag,np.matmul(U6dag,U5)))).real
						
						if(Sold<Snew):
							if(userandnums[trial]<=math.exp(Sold-Snew)):
								config[x,y,t,mu]=U1new;

						if(Sold>Snew):
								config[x,y,t,mu]=U1new;						

		U=config[5,5,5,1]
		Udag[0][0]=np.conj(U[0][0])
		Udag[0][1]=np.conj(U[1][0])
		Udag[1][0]=np.conj(U[0][1])
		Udag[1][1]=np.conj(U[1][1])

		V11[runnum]=V11[runnum]+np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig1,Udag))))).real*np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig1,Udag))))).real;
		V13[runnum]=V13[runnum]+np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig3,Udag))))).real*np.trace((np.matmul(sig1,np.matmul(U,np.matmul(sig3,Udag))))).real;

	print(bvals[runnum],V11[runnum]/usefultrials,V13[runnum]/usefultrials)

V11=V11/usefultrials
V13=V13/usefultrials

plt.plot(bvals,V11, 'r-')
plt.plot(bvals,V13, 'b-')
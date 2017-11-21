#Model from Huissman and Weissing 1995 Am. Nat. modified to use DOC inputs to change average light climate to estimate GPP -- this version is the original, loads are handled as only making contact with the epilimnion with all loads essentially going to the mixed layer
#Stuart Jones & Patrick Kelly 20  November 2017

library(deSolve)
#make dailyPAR data - just use some random incident light data
I0<-300

lightAtten<-function(z,I0,kD){
	Iz=I0*exp(-kD*z)
	return(Iz)
	}

# function for simulationg lake primary productivity as a function of DOC and P supply, lake surface area, and residence time
huisman<-function(t,y,params){
	with(as.list(c(y,params)),{
		
		# using Morris paper relationship, but symbols from Huisman
		kD=kA*A+kDOC*DOC	-0.05	#m-1; based on Morris paper

		# from a published paper -> I'll refind the paper eventually
		zmix=10^(-0.515*log10(DOC)+0.115*log10(2*sqrt(SA/pi))+0.991)
		if(zmix>zmax){zmix=zmax}
		
		V=SA*1e6*zmax
		Qin=V/365

		Izmix=lightAtten(z=zmix,I0=I0,kD=kD)
		
		# biomass specific growth integrated across mixed layer
		prod=(pA/(kD*zmix))*log((hA+I0)/(hA+Izmix))*(P/(P+mA))	# d-1
		
		dA.dt=A*prod-lA*A-v/zmix*A-Qin/(zmix*SA*1e6)*A	# mg C m-3
		dP.dt=Qin/(zmix*SA*1e6)*(Pin-P)+cA*lA*A*rec-cA*A*prod		# mg P m-3 (in epi)
		dDOC.dt=(Qin/(zmix*SA*1e6))*(DOCin-DOC)-decay*DOC				# g C m-3 (in epi)
		
		return(list(c(dA.dt,dP.dt,dDOC.dt)))
	})
}

## parameters
#SA=10		# km
#zmax=10			# m
#kDOC=0.22		# (m-1)/(g C m-3); from Morris paper; or 0.1-5.6 from Jager & Diehl
#kA=0.07/75/1000	# (m-1)/(mg C m-3); from Morris paper or 0.0003 from Jager & Diehl
#lA=0.1			# day-1
#pA=1		# day-1
#hA=100			# uE
#mA=3				# mg P m-3
#decay=0.001		# day-1
#Qin=SA*1e6*zmax/365	# m3
#Pin=5			# mg P m-3
#DOCin=5			# g C m-3
#cA=0.015			# (mg P) (mg C)-1; used to be 0.026
#V=SA*1e6*zmax	# m3
#v=0.1			# m d-1; sinking loss of algae
#rec=0.9			# proportion of P recycled from dead phytoplankton

times=1:10000

#SAs<-c(0.01,0.1,1,10)
SAs=c(0.01,0.1,1,10)
#SAs=0.01
DOCs<-seq(1,40,length.out=30)
#DOCs<-seq(1,40,length.out=12)
#CPs<-c(0.01,0.05,0.1)
CPs<-c(0.2,0.6,1.2)				# (g C m-3)*(mg P m-3)-1; rather than running DOC and P loads independently running three trajectories that differ in their load stoich
zmax=10

Ps<-seq(5,150,length.out=30)

# storing equilibrium state variables across simulations
storeAs=array(NA,dim=c(length(DOCs),length(CPs),length(SAs)))
storePs=storeAs
storeDOCs=storeAs

# loop through different surface areas, DOC loads, and load stoichs
for(k in 1:length(SAs)){
	for(j in 1:length(DOCs)){
		for(i in 1:length(CPs)){

#submittd
#parms=c(SA=SAs[k],zmax=10,kDOC=0.22,kA=0.00015,lA=0.1,pA=1,hA=35,mA=3,decay=0.001,Qin=SAs[k]*1e6*zmax/365,Pin=DOCs[j]/CPs[i],DOCin=DOCs[j],cA=0.0045,v=0.1,rec=0.9)
# used in jager & diehl
#parms=c(SA=SAs[k],zmax=10,kDOC=0.22,kA=0.00015,lA=0.1,pA=1,hA=100,mA=3,decay=0.001,Qin=SAs[k]*1e6*zmax/365,Pin=DOCs[j]/CPs[i],DOCin=DOCs[j],cA=0.015,v=0.1,rec=0.9)
# tweaking to get range of observations
parms=c(SA=SAs[k],zmax=10,kDOC=0.42,kA=0.00022,lA=0.1,pA=1.2,hA=55,mA=2,decay=0.001,Pin=DOCs[j]/CPs[i],DOCin=DOCs[j],cA=0.015,v=0.05,rec=0.99)


			# starting state variables
			n<-c(A=100,P=(DOCs[j]/CPs[i]),DOC=DOCs[j])
			
			# simulate with ode
			run=ode(y=n,times=times,func=huisman,parms=parms)
			
			# store equilibrium values
			storeAs[j,i,k]<-run[nrow(run),2]
			storePs[j,i,k]<-run[nrow(run),3]
			storeDOCs[j,i,k]<-run[nrow(run),4]
			
			print(c(k,i,j))
		}
	}
}

# calculate other simulation characteristics from stored equilibrium state variables
store_kD=parms[4]*storeAs+parms[3]*storeDOCs-0.05
store_zmix<-store_kD
for(i in 1:dim(store_zmix)[1]){
	for(j in 1:dim(store_zmix)[2]){
		for(k in 1:dim(store_zmix)[3]){
			store_zmix[i,j,k]<-10^(-0.515*log10(storeDOCs[i,j,k])+0.115*log10(2*sqrt(SAs[k]/pi))+0.991)

		}
	}
}
store_zmix[store_zmix>zmax]=zmax
storeTP=storeAs*parms[13]+storePs
storePP=store_kD*0
light.limit.d<-store_kD*0
nutrient.limit.d<-store_kD*0
storer<-store_kD*0
for(i in 1:dim(store_kD)[1]){
	for(j in 1:dim(store_kD)[2]){
		for(k in 1:dim(store_kD)[3]){	
				
		cP=storePs[i,j,k]
		
		kD=store_kD[i,j,k]

		zmix=store_zmix[i,j,k]
				
		Izmix=lightAtten(z=zmix,I0=I0,kD=kD)
		
		pA=parms[6]
		hA=parms[7]
		mA=parms[8]
		
		storer[i,j,k]<-(pA/(kD*zmix))*log((hA+I0)/(hA+Izmix))*(cP/(cP+mA))
		storePP[i,j,k]<-storeAs[i,j,k]*(pA/(kD*zmix))*log((hA+I0)/(hA+Izmix))*(cP/(cP+mA))	
		light.limit.d[i,j,k]=log((hA+I0)/(hA+Izmix))
		nutrient.limit.d[i,j,k]<-(cP/(cP+mA))
	
		}
	}
}

store_arealPP=storePP*store_zmix


# plotting
quartz()
plot(storeDOCs[,1,1],store_arealPP[,1,1],type='l',ylim=c(0,max(store_arealPP)),xlab="DOC (g C m-3)",ylab="Areal GPP (mg C m-2 day-1)",xlim=c(0,40))
lines(storeDOCs[,2,1],store_arealPP[,2,1],lty=2)
lines(storeDOCs[,3,1],store_arealPP[,3,1],lty=3)
lines(storeDOCs[,1,2],store_arealPP[,1,2],col='red')
lines(storeDOCs[,2,2],store_arealPP[,2,2],lty=2,col='red')
lines(storeDOCs[,2,3],store_arealPP[,3,2],lty=3,col='red')
legend('topleft',c('C:P=0.2, SA=0.01 km2','C:P=0.6, SA=0.01 km2','C:P=1.2, SA=0.01 km2', 'C:P=0.2, SA=10 km2','C:P=0.6, SA=10 km2','C:P=1.2, SA=10 km2'),lty=c(1:3,1:3),col=rep(c('black','red'),each=3),box.lty=0)

#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#define NUMBER_OF_LAYERS 4
#define NUMBER_OF_TRAINING_SET 4
double features[4][3]={
	{1.0,0.0,0.0},
	{1.0,0.0,1.0},
	{1.0,1.0,0.0},
	{1.0,1.0,1.0},
};
double labels[4][2]={
	{0.0,0.0},
	{1.0,0.0},
	{1.0,0.0},
	{0.0,1.0}
};
double randn (double mu, double sigma)
{
double U1, U2, W, mult;
static double X1, X2;
static int call = 0;
if (call == 1)
{
call = !call;
return fabs(mu + sigma * (double) X2);
}
do
{
U1 = -1 + ((double) rand () / RAND_MAX) * 2;
U2 = -1 + ((double) rand () / RAND_MAX) * 2;
W = pow (U1, 2) + pow (U2, 2);
}while (W >= 1 || W == 0);
mult = sqrt ((-2 * log (W)) / W);
X1 = U1 * mult;
X2 = U2 * mult;
call = !call;
return fabs(mu + sigma * (double) X1);
}
int main()
{
int* layers;
double* neurons;
double* dendrites;
double* gradients;
double** lb;
double** db;
double** gb;
int numberOfNeurons,numberOfDendrites;
int i,j,k,l,m,n,offset,cycle,dataSetExNo,layerNum;
double *ref;
double mu,sigma;
double *currentNeuron,*currentDendron,*currentGradient,*prevLayer;
double sum;
double error,acc;
double deltaWeight;
double* tmp;
layers=(int*)malloc(sizeof(int)*NUMBER_OF_LAYERS);
layers[0]=2+1;
layers[1]=4+1;
layers[2]=3+1;
layers[3]=2;
for(i=0,numberOfNeurons=0;i<NUMBER_OF_LAYERS;++i)
{
numberOfNeurons+=layers[i];
}
neurons=(double*)malloc(sizeof(double*)*numberOfNeurons);
gradients=(double*)malloc(sizeof(double*)*numberOfNeurons);
for(i=0,numberOfDendrites;i<NUMBER_OF_LAYERS;++i)
{
numberOfDendrites+=(layers[i]*layers[i+1]);
}
dendrites=(double*)malloc(sizeof(double)*numberOfDendrites);
lb=(double**)malloc(sizeof(double*)*NUMBER_OF_LAYERS);
for(i=0,offset=0;i<NUMBER_OF_LAYERS;++i)
{
lb[i]=neurons+offset;
offset+=layers[i];
}
db=(double**)malloc(sizeof(double*)*(NUMBER_OF_LAYERS-1));
for(i=0,offset=0;i<NUMBER_OF_LAYERS;++i)
{
db[i]=dendrites+offset;
offset+=(layers[i]*layers[i+1]);
}
gb=(double**)malloc(sizeof(double*)*(NUMBER_OF_LAYERS));
for(i=0,offset=0;i<NUMBER_OF_LAYERS;++i)
{
gb[i]=gradients+offset;
offset+=layers[i];
}
mu=0.0;
sigma=0.1;
double val=0.0;
for(i=0,tmp=dendrites;i<numberOfDendrites;++i,++val)
{
*tmp=randn(mu,sigma);
++tmp;
}


cycle=0;

while(1)
{
if(cycle%10000==0)printf("acc -->>>>>>>>>%f\n",acc);
acc=0.0;
for(dataSetExNo=0;dataSetExNo<NUMBER_OF_TRAINING_SET;++dataSetExNo)
{

for(i=0;i<layers[0];++i)
{
lb[0][i]=features[dataSetExNo][i];
}
//Feed-Forward starts
for(layerNum=1;layerNum<=(NUMBER_OF_LAYERS-1);++layerNum)
{
currentNeuron=lb[layerNum];
currentDendron=db[layerNum-1];
m=layers[layerNum];
for(i=0;i<m;++i)
{
prevLayer=lb[layerNum-1];
n=layers[layerNum-1];
sum=0.0;
for(j=0;j<n;++j)
{
sum+=((*prevLayer)*(*currentDendron));
++currentDendron;
++prevLayer;
}
*currentNeuron=1/(1+(exp((-1.0)*sum)));
++currentNeuron;
}
if(layerNum!=NUMBER_OF_LAYERS-1)*lb[layerNum]=1.0;
}

if(cycle%1000000==0)
{
printf("Input -> %f %f\n",features[dataSetExNo][1],features[dataSetExNo][2]);
printf("Output ->%f %f\n",lb[NUMBER_OF_LAYERS-1][0],lb[NUMBER_OF_LAYERS-1][1]);
printf("Target ->%f %f\n",labels[dataSetExNo][0],labels[dataSetExNo][1]);
}

//Feed-Forward ends

//Calculating output layer gradients
sum=0.0;
currentNeuron=lb[NUMBER_OF_LAYERS-1];
tmp=gb[NUMBER_OF_LAYERS-1];
m=layers[NUMBER_OF_LAYERS-1];
for(i=0;i<m;++i)
{
error=labels[dataSetExNo][i]-(*currentNeuron);
*tmp=error*(*(currentNeuron))*(1.0-(*currentNeuron));
sum+=(error*error);
++tmp;
++currentNeuron;
}
acc+=sqrt(sum/layers[NUMBER_OF_LAYERS-1]);

//Calculating hidden layer gradients
for(layerNum=NUMBER_OF_LAYERS-2;layerNum>0;--layerNum)
{
m=layers[layerNum];
currentGradient=gb[layerNum];
currentDendron=db[layerNum];
currentNeuron=lb[layerNum];
for(i=0;i<m;++i)
{
n=layers[layerNum+1];
tmp=gb[layerNum+1];
sum=0.0;
for(j=0;j<n;++j)
{
sum+=(*tmp)*(*currentDendron);
++tmp;
++currentDendron;
}
*currentGradient=sum*((*currentNeuron)*(1.0-(*currentNeuron)));
++currentGradient;
++currentNeuron;
}
}


for(layerNum=1;layerNum<NUMBER_OF_LAYERS;++layerNum)
{
currentGradient=gb[layerNum];
currentDendron=db[layerNum-1];
m=layers[layerNum];
for(i=0;i<m;++i)
{
n=layers[layerNum-1];
currentNeuron=lb[layerNum-1];
for(j=0;j<n;++j)
{
deltaWeight=(0.01)*(*currentGradient)*(*currentNeuron);
*currentDendron+=deltaWeight;
++currentNeuron;
++currentDendron;
}
++currentGradient;
}

}


}//Biggest for loop ever
cycle++;
}//Biggest while loop ever 





free(neurons);
free(gradients);
free(dendrites);
free(lb);
free(gb);
free(db);
free(layers);

return 0;
}

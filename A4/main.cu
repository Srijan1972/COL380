#include<bits/stdc++.h>
using namespace std;
#define MAX_VAL 4294967295LL

__global__ void matrixMul(uint16_t* matA,uint16_t* matB,uint32_t* matC,bool* existA,bool* existB,bool* existC,int n,int m,int tot){
    int elem = blockIdx.x*blockDim.x + threadIdx.x;
    if(elem<n*n){
        int x = elem/n;
        int y = elem%n;
        int i = x/m;
        int j = y/m;
        int a = x%m;
        int b = y%m;
        long long temp = 0;
        for(int k=0;k<tot;k++){
            bool A = existA[i*tot+k];
            bool B = existB[k*tot+j];
            if(A && B){
                int baseA = i*m*n + k*m*m;
                int baseB = k*m*n + j*m*m;
                for(int c=0;c<m;c++){
                    long long t1 = matA[baseA+a*m+c];
                    long long t2 = matB[baseB+c*m+b];
                    temp += t1*t2;
                }
            }
        }
        int baseC = i*m*n + j*m*m;
        long long t3 = min(MAX_VAL,temp);
        uint32_t tempC = t3;
        matC[baseC+a*m+b] = tempC;
        if(tempC!=0) existC[i*tot+j] = 1;
    }
}

int main(int argc,char **argv){
    if(argc!=4){
        cout<<"Invalid format: expected exec input1 input2 output\n";
        return 1;
    }
    int n,m,kA,kB,kC;
    FILE* inpA = fopen(argv[1],"rb");
    n = getw(inpA);
    m = getw(inpA);
    int tot = n/m;
    kA = getw(inpA);
    uint16_t *host_matA = new uint16_t[n*n]();
    uint16_t *host_matB = new uint16_t[n*n]();
    uint16_t *matA,*matB;
    cudaMalloc(&matA,n*n*sizeof(uint16_t));
    cudaMalloc(&matB,n*n*sizeof(uint16_t));
    uint32_t *host_matC = new uint32_t[n*n]();
    uint32_t *matC;
    cudaMalloc(&matC,n*n*sizeof(uint32_t));
    bool *host_existA = new bool[tot*tot]();
    bool *host_existB = new bool[tot*tot]();
    bool *host_existC = new bool[tot*tot]();
    bool *existA,*existB,*existC;
    cudaMalloc(&existA,tot*tot*sizeof(bool));
    cudaMalloc(&existB,tot*tot*sizeof(bool));
    cudaMalloc(&existC,tot*tot*sizeof(bool));
    for(int it=0;it<kA;it++){
        int i = getw(inpA);
        int j = getw(inpA);
        host_existA[i*tot+j] = 1;
        int t = i*m*n+j*m*m;
        for(int u = 0;u<m*m;u++){
            uint8_t low = getc(inpA);
            uint8_t upp = getc(inpA);
            host_matA[t+u] = (upp << 8) + low;
        }
    }
    cudaMemcpy(existA,host_existA,tot*tot*sizeof(bool),cudaMemcpyHostToDevice);
    cudaMemcpy(matA,host_matA,n*n*sizeof(uint16_t),cudaMemcpyHostToDevice);
    fclose(inpA);
    FILE* inpB = fopen(argv[2],"rb");
    assert(n==getw(inpB));
    assert(m==getw(inpB));
    kB = getw(inpB);
    for(int it=0;it<kB;it++){
        int i = getw(inpB);
        int j = getw(inpB);
        host_existB[i*tot+j] = 1;
        int t = i*m*n+j*m*m;
        for(int u = 0;u<m*m;u++){
            uint8_t low = getc(inpB);
            uint8_t upp = getc(inpB);
            host_matB[t+u] = (upp << 8) + low;
        }
    }
    cudaMemcpy(existB,host_existB,tot*tot*sizeof(bool),cudaMemcpyHostToDevice);
    cudaMemcpy(matB,host_matB,n*n*sizeof(uint16_t),cudaMemcpyHostToDevice);
    fclose(inpB);
    matrixMul<<<(n*n+1023)/1024,1024>>>(matA,matB,matC,existA,existB,existC,n,m,tot);
    cudaDeviceSynchronize();
    cudaMemcpy(host_existC,existC,tot*tot*sizeof(bool),cudaMemcpyDeviceToHost);
    cudaMemcpy(host_matC,matC,n*n*sizeof(uint32_t),cudaMemcpyDeviceToHost);
    kC = 0;
    for(int i=0;i<tot*tot;i++){
        if(host_existC[i]) kC++;
    }
    FILE* outC = fopen(argv[3],"wb");
    putw(n,outC);
    putw(m,outC);
    putw(kC,outC);
    for(int w=0;w<tot*tot;w++){
        if(host_existC[w]){
            int i = w/tot;
            int j = w%tot;
            putw(i,outC);
            putw(j,outC);
            int base = i*m*n + j*m*m;
            for(int u=0;u<m*m;u++) putw(host_matC[base+u],outC);
        }
    }
    fclose(outC);
    cudaFree(matA);
    cudaFree(existA);
    cudaFree(matB);
    cudaFree(existB);
    cudaFree(matC);
    cudaFree(existC);
    return 0;
}

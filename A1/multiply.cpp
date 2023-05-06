#include "multiply.hpp"
#include "library.hpp"
#include <omp.h>
#define MAX_VAL 65535

// * Perform cumulative "multiplacation" of two blocks
void basic_mult(uint8_t* m1,uint8_t* m2,int* res,int m){
    for(int i=0;i<m;i++){
        for(int k=0;k<m;k++){
            for(int j=0;j<m;j++){
                res[i*m+j] = Outer(res[i*m+j],Inner(m1[i*m+k],m2[k*m+j]));
            }
        }
    }
}

// This function reads from the input file, computes square and writes it to the output file
void matrix_square(string infile,string outfile)
{
    // * Reading Input
    ifstream inp;
    inp.open(infile);
    int n;
    inp.read(reinterpret_cast<char*>(&n), sizeof(n));
    int m;
    inp.read(reinterpret_cast<char*>(&m), sizeof(m));
    int k;
    inp.read(reinterpret_cast<char*>(&k), sizeof(k));
    map<int,map<int,uint8_t*>> inp_mat;
    for(int l=0;l<k;l++){
        int i,j;
        inp.read(reinterpret_cast<char*>(&i), sizeof(i));
        inp.read(reinterpret_cast<char*>(&j), sizeof(j));
        inp_mat[i][j] = new uint8_t[m*m];
        if(i!=j) inp_mat[j][i] = new uint8_t[m*m];
        for(int a=0;a<m;a++){
            for(int b=0;b<m;b++){
                uint8_t c;
                inp.read(reinterpret_cast<char*>(&c), sizeof(c));
                inp_mat[i][j][a*m+b] = c;
                if(i!=j) inp_mat[j][i][b*m+a] = c;
            }
        }
    }
    // * Computing Matrix
    int tot=n/m;
    map<int,map<int,int*>> out_mat;
    #pragma omp parallel
    {
        #pragma omp single
        for(int i=0;i<tot;i++){
            for(int l=0;l<tot;l++){
                uint8_t *arg1;
                auto it1 = inp_mat.find(i);
                if(it1==inp_mat.end()) arg1 = NULL;
                else{
                        auto f = it1->second.find(l);
                        if(f!=it1->second.end()) arg1 = f->second;
                        else arg1 = NULL;
                }
                if(arg1!=NULL){
                    #pragma omp task
                    {
                        for(int j=i;j<tot;j++){
                            uint8_t *arg2;
                            auto it2 = inp_mat.find(l);
                            if(it2==inp_mat.end()) arg2 = NULL;
                            else{
                                auto f = it2->second.find(j);
                                if(f!=it2->second.end()) arg2 = f->second;
                                else arg2 = NULL;
                            }
                            if(arg2){
                                int* ou = new int[m*m];
                                for(int a=0;a<m*m;a++) ou[a]=0;
                                basic_mult(arg1,arg2,ou,m);
                                #pragma omp critical
                                {
                                    if(out_mat[i][j]==NULL){
                                        out_mat[i][j] = ou;
                                    }
                                    else{
                                        for(int a=0;a<m*m;a++){
                                            out_mat[i][j][a] = Outer(out_mat[i][j][a],ou[a]);
                                        }
                                    }
                                }
                                // #pragma omp critical
                            }
                        }
                    }
                    // #pragma omp taskwait
                }
            }
        }
        #pragma omp taskwait
    }
    // * Writing Output
    auto it = out_mat.begin();
    int K=0;
    while(it!=out_mat.end()){
        K+=(*it).second.size();
        it++;
    }

    ofstream out;
    out.open(outfile,ios::binary);
    out.write(reinterpret_cast<char*>(&n), sizeof(n));
    out.write(reinterpret_cast<char*>(&m), sizeof(m));
    out.write(reinterpret_cast<char*>(&K), sizeof(K));
    for(int i=0;i<tot;i++){
        for(int j=i;j<tot;j++){
            if(out_mat[i][j]){
                out.write(reinterpret_cast<char*>(&i), sizeof(i));
                out.write(reinterpret_cast<char*>(&j), sizeof(j));
                uint16_t c;
                for(int a=0;a<m*m;a++){
                    int temp = out_mat[i][j][a];
                    temp = min(temp,MAX_VAL);
                    c = temp;
                    out.write(reinterpret_cast<char*>(&c), sizeof(c));
                }
            }
        }
    }
}

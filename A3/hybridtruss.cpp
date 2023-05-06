#include<bits/stdc++.h>
#include<omp.h>
#include<mpi.h>
using namespace std;

bool compareVert(vector<int> &numAdj,int i,int j){
    int t1 = numAdj[i];
    int t2 = numAdj[j];
    return (t1<t2 || (t1==t2 && i<j));
}

bool update(pair<int,int> e,int vold,int vnew,map<pair<int,int>,int>& T,map<pair<int,int>,int>& g,map<pair<int,int>,map<int,int>>& h){
    int temp = T[e];
    bool flag = false;
    if(vold>=temp && vnew<temp){
        g[e]--;
        h[e][vnew]++;
    }
    else if(vold<temp && vnew<temp){
        h[e][vnew]++;
        h[e][vold]--;
    }
    if(g[e]<temp-2){
        temp--;
        g[e]+=h[e][temp];
        flag = true;
    }
    T[e] = temp;
    return flag;
} 

void BFS(string outputpath,int n,vector<set<int>>& graph){
    FILE* f = fopen(outputpath.c_str(),"a");
    vector<bool> visited(n,false);
    vector<set<int>> truss;
    for(int v=0;v<n;v++){
        if(visited[v] || graph[v].size()==0) continue;
        queue<int> vert;
        vert.push(v);
        set<int> t;
        while(!vert.empty()){
            int u = vert.front();
            vert.pop();
            if(visited[u] || graph[u].size()==0) continue;
            visited[u] = true;
            t.insert(u);
            for(auto w:graph[u]) vert.push(w);
        }
        truss.push_back(t);
    }
    fprintf(f,"%d\n",truss.size()>0);
    if(truss.size()){
        fprintf(f,"%ld\n",truss.size());
        for(auto t:truss){
            for(int v:t) fprintf(f,"%d ",v);
            fprintf(f,"\n");
        }
    }
    fclose(f);
}

class InputParser{
    private:
        vector<string> args;
    public:
        InputParser(int argc,char* argv[]){
            args.clear();
            for(int i=1;i<argc;i++) args.push_back(string(argv[i]));
        }
        bool findOption(string arg){
            auto it = args.begin();
            while(it!=args.end()){
                string opt = *it;
                if(opt.length()>arg.length() && opt.substr(0,arg.length())==arg) return true;
                it++;
            }
            return false;
        }
        string getOption(string arg){
            auto it = args.begin();
            while(it!=args.end()){
                string opt = *it;
                if(opt.length()>arg.length() && opt.substr(0,arg.length())==arg) return opt.substr(arg.length()+1);
                it++;
            }
            return "";
        }
};

int main(int argc,char* argv[]){
    int why;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&why);
    assert(MPI_THREAD_MULTIPLE==why);
    InputParser* cmd = new InputParser(argc,argv);
    int taskid=1;
    if(cmd->findOption("--taskid")){
        int temp = stoi(cmd->getOption("--taskid"));
        assert(temp>=1 && temp<=2);
        taskid = temp;
    }
    string inputpath,headerpath,outputpath;
    bool verbose=0;
    int startk=-1,endk=-1,p=-1;
    if(!cmd->findOption("--inputpath")){
        cerr<<"Input path parameter not passed\n";
        return 1;
    }
    else inputpath = cmd->getOption("--inputpath");
    if(!cmd->findOption("--headerpath")){
        cerr<<"Header path parameter not passed\n";
        return 1;
    }
    else headerpath = cmd->getOption("--headerpath");
    if(!cmd->findOption("--outputpath")){
        cerr<<"Output path parameter not passed\n";
        return 1;
    }
    else outputpath = cmd->getOption("--outputpath");
    if(!cmd->findOption("--verbose")){
        cerr<<"verbose parameter not passed\n";
        return 1;
    }
    else{
        int temp = stoi(cmd->getOption("--verbose"));
        assert(temp>=0 && temp<=1);
        verbose = temp;
    }
    if(cmd->findOption("--startk")) startk = stoi(cmd->getOption("--startk"));
    if(cmd->findOption("--endk")) endk = stoi(cmd->getOption("--endk"));
    if(cmd->findOption("--p")) p = stoi(cmd->getOption("--p"));
    if(taskid==1 && startk<0){
        cerr<<"startk required for task 1\n";
        return 1;
    }
    if(endk<0){
        cerr<<"endk parameter not passed\n";
        return 1;
    }
    if(taskid==2 && p<0){
        cerr<<"p required for task 2\n";
        return 1;
    }
    int numProcs,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status status;
    int n,m;
    FILE* inp = fopen(inputpath.c_str(),"rb");
    FILE* head = fopen(headerpath.c_str(),"rb");
    n = getw(inp);
    m = getw(inp);
    vector<int> offsets;
    vector<set<int>> adjList;
    vector<int> deg;
    for(int i=0;i<n;i++){
        offsets.push_back(getw(head));
        fseek(inp,offsets[i]+4,SEEK_SET);
        deg.push_back(getw(inp));
    }
    fclose(head);
    for(int i=0;i<n;i++){
        fseek(inp,offsets[i]+8,SEEK_SET);
        set<int> temp;
        for(int j=0;j<deg[i];j++){
            int v = getw(inp);
            if(compareVert(deg,i,v)) temp.insert(v);
        }
        adjList.push_back(temp);
    }
    fclose(inp);
    // map<pair<int,int>,int> supp;
    map<pair<int,int>,map<int,int>> tris;
    vector<vector<int>> triSend(numProcs,vector<int>());
    #pragma omp parallel for
    for(int u=rank;u<n;u+=numProcs){
        #pragma omp task
        {
            for(int v:adjList[u]){
                int proc = v%numProcs;
                vector<int> comm;
                set_intersection(adjList[u].begin(),adjList[u].end(),adjList[v].begin(),adjList[v].end(),inserter(comm,comm.begin()));
                if(proc==rank){
                    for(int w:comm){
                        #pragma omp critical (suppset)
                        {
                            // supp[{u,v}]++;
                            // supp[{v,w}]++;
                            // supp[{u,w}]++;
                            tris[{u,v}][w] = INT_MAX;
                            tris[{v,w}][u] = INT_MAX;
                            tris[{u,w}][v] = INT_MAX;
                        }
                    }
                }
                else{
                    for(int w:comm){
                        #pragma omp critical (suppset)
                        {
                            // supp[{u,v}]++;
                            // supp[{u,w}]++;
                            tris[{u,v}][w] = INT_MAX;
                            tris[{u,w}][v] = INT_MAX;
                       }
                        #pragma omp critical (suppsend)
                        {
                            triSend[proc].push_back(v);
                            triSend[proc].push_back(w);
                            triSend[proc].push_back(u);
                        }
                    }
                }
            }
        }
    }
    int* sizes = new int[numProcs];
    for(int i=0;i<numProcs;i++) sizes[i] = triSend[i].size();
    int* recs = new int[numProcs];
    MPI_Alltoall(sizes,1,MPI_INT,recs,1,MPI_INT,MPI_COMM_WORLD);
    vector<int> dat;
    for(auto it:triSend){
        for(int a:it) dat.push_back(a);
    }
    int* off = new int[numProcs];
    off[0] = 0;
    for(int i=1;i<numProcs;i++) off[i]=off[i-1]+sizes[i-1];
    int* roff = new int[numProcs];
    roff[0] = 0;
    for(int i=1;i<numProcs;i++) roff[i]=roff[i-1]+recs[i-1];
    int* triRecv = new int[roff[numProcs-1]+recs[numProcs-1]];
    MPI_Alltoallv(dat.data(),sizes,off,MPI_INT,triRecv,recs,roff,MPI_INT,MPI_COMM_WORLD);
    #pragma omp parallel for
    for(int l=0;l<roff[numProcs-1]+recs[numProcs-1];l+=3){
        #pragma omp task
        {
            int v = triRecv[l];
            int w = triRecv[l+1];
            int u = triRecv[l+2];
            #pragma omp critical
            {
                // supp[{v,w}]++;
                tris[{v,w}][u] = INT_MAX;
            }
        }
    }
    map<pair<int,int>,int> T;
    map<pair<int,int>,int> g;
    map<pair<int,int>,map<int,int>> h;
    for(int u=rank;u<n;u+=numProcs){
        for(int v:adjList[u]){
            int supp = tris[{u,v}].size();
            T[{u,v}] = supp+2;
            g[{u,v}] = supp;
            for(int j=0;j<supp+2;j++) h[{u,v}][j] = 0;
        }
    }
    map<pair<int,int>,int> gam = T;
    int kmin=INT_MAX,kmax=INT_MIN;
    for(auto p:T){
        kmin=min(kmin,p.second);
        kmax=max(kmax,p.second);
    }
    int k=kmin;
    set<pair<int,int>> Act;
    int progress = 1;
    for(auto p:T){
        if(p.second==k) Act.insert(p.first);
    }
    while(Act.size() || k<kmax || progress){
        set<pair<int,int>> changed;
        int done = 0,flag = 0;
        for(auto e:Act){
            int u = e.first;
            int v = e.second;
            int su = T[e];
            for(auto kv:tris[e]){
                int w = kv.first;
                int vold = kv.second;
                if(su<vold){
                    int a=u,b=w,c=v,d=w;
                    if(compareVert(deg,b,a)) swap(a,b);
                    if(compareVert(deg,d,c)) swap(c,d);
                    int vnew = su;
                    tris[{u,v}][w] = vnew;
                    if(a%numProcs==rank){
                        tris[{a,b}][v] = vnew;
                        if(update({a,b},vold,vnew,T,g,h) && T[{a,b}]>=kmin && T[{a,b}]<=k)changed.insert({a,b});
                    }
                    else{
                        int trip[4] = {a,b,v,vnew};
                        MPI_Send(trip,4,MPI_INT,a%numProcs,0,MPI_COMM_WORLD);
                    }
                    if(c%numProcs==rank){
                        tris[{c,d}][u] = vnew;
                        if(update({c,d},vold,vnew,T,g,h) && T[{c,d}]>=kmin && T[{c,d}]<=k) changed.insert({c,d});
                    }
                    else{
                        int trip[4] = {c,d,u,vnew};
                        MPI_Send(trip,4,MPI_INT,c%numProcs,0,MPI_COMM_WORLD);
                    }
                }
            }
        }
        int fin[4] = {-1,-1,-1,-1};
        for(int i=0;i<numProcs;i++){
            if(i!=rank) MPI_Send(fin,4,MPI_INT,i,0,MPI_COMM_WORLD);
        }
        done++;
        while(done<numProcs){
            MPI_Iprobe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&flag,&status);
            if(flag){
                int rec[4];
                MPI_Recv(rec,4,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD,&status);
                if(rec[0]<0) done++;
                else{
                    int u = rec[0],v = rec[1],w = rec[2];
                    int vnew = rec[3];
                    int vold = tris[{u,v}][w];
                    tris[{u,v}][w] = vnew;
                    if(update({u,v},vold,vnew,T,g,h) && T[{u,v}]>=kmin && T[{u,v}]<=k) changed.insert({u,v});
                }
            }
        }
        Act = changed;
        if(Act.size()==0 && k==kmax) progress = 0;
        int temp;
        MPI_Allreduce(&progress,&temp,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        progress = temp;
        if(k!=kmax){
            k++;
            for(auto p:T){
                if(p.second==k) Act.insert(p.first);
            }
        }
    }
    int mx = 0;
    for(auto p:T) mx = max(mx,p.second);
    int mtruss;
    MPI_Allreduce(&mx,&mtruss,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(taskid==1){
        if(!verbose){
            if(rank==0){
                FILE* out = fopen(outputpath.c_str(),"w");
                for(int l=startk;l<=endk;l++){
                    if(l+2<=mtruss) fprintf(out,"1 ");
                    else fprintf(out,"0 ");
                }
                fprintf(out,"\n");
                fclose(out);
            }
        }
        else{
        int flag = 0;
        int edge[2];
        int k = startk;
        int req = startk + 2;
        int fin = endk + 2;
        if(rank==0){
            FILE* out = fopen(outputpath.c_str(),"w");
            fclose(out);
            vector<set<int>> graph(n,set<int>());
            for(auto it=T.begin();it!=T.end();it++){
                graph[it->first.first].insert(it->first.second);
                graph[it->first.second].insert(it->first.first);
            }
            int done = 1;
            while(done<numProcs){
                MPI_Iprobe(MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&flag,&status);
                if(flag){
                    MPI_Recv(edge,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD,&status);
                    if(edge[0]<0) done++;
                    else{
                        graph[edge[0]].insert(edge[1]);
                        graph[edge[1]].insert(edge[0]);
                    }
                }
            }
            while(req<=fin){
                for(int i=1;i<numProcs;i++) MPI_Send(&req,1,MPI_INT,i,2,MPI_COMM_WORLD);
                done = 1;
                while(done<numProcs){
                    MPI_Iprobe(MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&flag,&status);
                    if(flag){
                        MPI_Recv(edge,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD,&status);
                        if(edge[0]<0) done++;
                        else{
                            graph[edge[0]].erase(edge[1]);
                            graph[edge[1]].erase(edge[0]);
                        }
                    }
                }
                for(auto it=T.begin();it!=T.end();){
                    if(it->second<req){
                        graph[it->first.first].erase(it->first.second);
                        graph[it->first.second].erase(it->first.first);
                        T.erase(it++);
                    }
                    else ++it;
                }
                BFS(outputpath,n,graph);
                req++;
            }
            req = -1;
            for(int i=1;i<numProcs;i++) MPI_Send(&req,1,MPI_INT,i,2,MPI_COMM_WORLD);
        }
        else{
            for(auto it=T.begin();it!=T.end();it++){
                edge[0] = it->first.first; edge[1] = it->first.second;
                MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
            }
            edge[0] = -1; edge[1] = -1;
            MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
            while(1){
                MPI_Iprobe(0,2,MPI_COMM_WORLD,&flag,&status);
                if(flag){
                    MPI_Recv(&req,1,MPI_INT,0,2,MPI_COMM_WORLD,&status);
                    if(req<0) break;
                    for(auto it=T.begin();it!=T.end();){
                        if(it->second<req){
                            edge[0] = it->first.first; edge[1] = it->first.second;
                            MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
                            T.erase(it++);
                        }
                        else ++it;
                    }
                    edge[0] = -1; edge[1] = -1;
                    MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
                }
            }
        }
    }
    }
    if(taskid==2){
        int req = endk + 2;
        int edge[2];
        if(rank==0){
            vector<set<int>> graph(n,set<int>());
            vector<int> pVals(n,0);
            map<int,vector<int>> pMap;
            if(verbose){
                for(int i=0;i<n;i++){
                    vector<int> v;
                    pMap[i] = v;
                }
            }
            for(auto it=T.begin();it!=T.end();it++){
                if(it->second>=req){
                    graph[it->first.first].insert(it->first.second);
                    graph[it->first.second].insert(it->first.first);
                }
            }
            int done = 1,flag;
            while(done<numProcs){
                MPI_Iprobe(MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&flag,&status);
                if(flag){
                    MPI_Recv(edge,2,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD,&status);
                    if(edge[0]<0) done++;
                    else{
                        graph[edge[0]].insert(edge[1]);
                        graph[edge[1]].insert(edge[0]);
                    }
                }
            }
            vector<set<int>> truss;
            vector<bool> visited(n,false);
            for(int i=0;i<n;i++){
                if(visited[i] || graph[i].size()==0) continue;
                queue<int> it;
                it.push(i);
                set<int> t;
                while(!it.empty()){
                    int u = it.front();
                    it.pop();
                    if(visited[u] || graph[u].size()==0) continue;
                    visited[u] = true;
                    t.insert(u);
                    for(auto v:graph[u]) it.push(v);
                }
                truss.push_back(t);
            }
            for(int i=0;i<n;i++){
                for(auto v:adjList[i]) adjList[v].insert(i);
            }
            int iter = 0;
            for(auto t:truss){
                set<int> upd;
                for(auto u:t){
                    for(auto v:adjList[u]) upd.insert(v);
                }
                for(int v:upd){
                    if(!verbose) pVals[v]++;
                    else pMap[v].push_back(iter);
                }
                iter++;
            }
            vector<int> influencers;
            if(!verbose){
                for(int i=0;i<n;i++){
                    if(pVals[i]>=p) influencers.push_back(i);
                }
                FILE* out = fopen(outputpath.c_str(),"w");
                fprintf(out,"%ld\n",influencers.size());
                for(int i=0;i<influencers.size();i++) fprintf(out,"%d ",influencers[i]);
                if(influencers.size()) fprintf(out,"\n");
                fclose(out);
            }
            else{
                for(int i=0;i<n;i++){
                    if(pMap[i].size()>=p) influencers.push_back(i);
                }
                FILE* out = fopen(outputpath.c_str(),"w");
                fprintf(out,"%ld\n",influencers.size());
                for(int i=0;i<influencers.size();i++){
                    fprintf(out,"%d\n",influencers[i]);
                    for(int it:pMap[influencers[i]]){
                        for(int u:truss[it]) fprintf(out,"%d ",u);
                    }
                    fprintf(out,"\n");
                }
                fclose(out);
            }
        }
        else{
            for(auto it=T.begin();it!=T.end();it++){
                if(it->second>=req){
                    edge[0] = it->first.first; edge[1] = it->first.second;
                    MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
                }
            }
            edge[0] = -1; edge[1] = -1;
            MPI_Send(edge,2,MPI_INT,0,1,MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}
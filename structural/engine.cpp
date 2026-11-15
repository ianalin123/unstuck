// Local search engine written for this project. No network or API calls.
// Prior art: CPro1 conflict/tabu search and ARW-style feasible MIS search.
// Search outputs are hypotheses; fixed_verifier.py is the independent authority.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <set>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace std;
using Clock=chrono::steady_clock;

struct Data {
    int n,N,M; vector<vector<int>> d,parents; vector<int> runs,weight,vt;
    Data(int length):n(length),N(1<<n),M(1<<(n-2)),d(N),parents(M),runs(N),weight(N),vt(N) {
        for(int w=0;w<N;++w) {
            auto &row=d[w];
            for(int i=0;i<n-1;++i) {
                int first=(w&((1<<i)-1))|((w>>(i+1))<<i);
                for(int j=i+1;j<n;++j) {
                    int k=j-1;row.push_back((first&((1<<k)-1))|((first>>(k+1))<<k));
                }
            }
            sort(row.begin(),row.end());row.erase(unique(row.begin(),row.end()),row.end());
            for(int c:row)parents[c].push_back(w);
            runs[w]=1;weight[w]=__builtin_popcount((unsigned)w);
            for(int k=0;k<n;++k) {
                vt[w]+=(k+1)*((w>>k)&1);
                if(k&&((w>>k)&1)!=((w>>(k-1))&1))++runs[w];
            }
            vt[w]%=n+1;
        }
        // Every binary (n-2)-word has the same number 1+n+C(n,2) of n-supersequences.
        for(auto &row:parents)assert((int)row.size()==1+n+n*(n-1)/2);
    }
};

struct Search {
    Data &g;mt19937_64 rng;vector<int> code,position,occ,tabu;
    int cost=0,bestCost=1000000,bestValid=0;long long iteration=0,evaluations=0;
    vector<int> penalties;
    vector<char> initialMask;set<pair<int,int>> archiveBins;
    Clock::time_point begin;double seconds;string out;vector<int> bestCode,bestState;
    Search(Data &data,int seed,double budget,string directory):g(data),rng(seed),position(g.N,-1),occ(g.M),tabu(g.N),penalties(g.M,1),begin(Clock::now()),seconds(budget),out(directory){}
    double elapsed(){return chrono::duration<double>(Clock::now()-begin).count();}
    bool timeLeft(){return elapsed()<seconds;}
    int random(int bound){return rng()%bound;}
    double unit(){return generate_canonical<double,53>(rng);}
    void add(int w){assert(position[w]<0);position[w]=code.size();code.push_back(w);for(int c:g.d[w]){if(occ[c]++)++cost;}}
    void eraseIndex(int i){int w=code[i];for(int c:g.d[w]){if(occ[c]-->1)--cost;}position[code.back()]=i;code[i]=code.back();code.pop_back();position[w]=-1;}
    int conflicts(int w){int k=0;for(int c:g.d[w])k+=occ[c]>0;return k;}
    int responsibility(int w){int k=0;for(int c:g.d[w])k+=occ[c]>1;return k;}
    void load(const vector<int>&words){fill(position.begin(),position.end(),-1);fill(occ.begin(),occ.end(),0);code.clear();cost=0;for(int w:words)add(w);}
    string binary(int w){string s(g.n,'0');for(int k=0;k<g.n;++k)s[g.n-1-k]='0'+((w>>k)&1);return s;}
    void save(const string&name,const vector<int>&words,int conflictCost){
        ofstream f(out+"/"+name);f<<"{\"n\":"<<g.n<<",\"conflict_cost\":"<<conflictCost<<",\"words\":[";
        for(size_t i=0;i<words.size();++i){if(i)f<<',';f<<'\"'<<binary(words[i])<<'\"';}f<<"]}\n";
    }
    void log(const string&event,const string&extra=""){
        cout<<"{\"event\":\""<<event<<"\",\"seconds\":"<<elapsed()<<",\"iteration\":"<<iteration<<",\"size\":"<<code.size()<<",\"cost\":"<<cost<<",\"best_valid\":"<<bestValid<<",\"best_cost\":"<<bestCost<<",\"evaluations\":"<<evaluations<<extra<<"}"<<endl;
    }
    bool record(){
        if(!initialMask.empty()&&cost<=20){
            int novel=0;for(int w:code)novel+=!initialMask[w];
            auto bin=make_pair(cost,novel/8);
            if(archiveBins.insert(bin).second)save("archive/state-c"+to_string(cost)+"-d"+to_string(novel/8)+".json",code,cost);
        }
        if(cost<bestCost){
            bestCost=cost;bestState=code;save("best-state.json",bestState,cost);log("lower_conflict_state");
            if(cost>0){
                auto state=code;
                while(cost>0){
                    int chosen=0,highest=-1;
                    for(int i=0;i<(int)code.size();++i){int r=responsibility(code[i]);if(r>highest){highest=r;chosen=i;}}
                    eraseIndex(chosen);
                }
                fillGreedy(0,16);
                if((int)code.size()>bestValid){bestValid=code.size();bestCode=code;save("best-valid.json",bestCode,0);log("feasible_extraction");}
                load(state);
            }
        }
        if(cost==0&&(int)code.size()>bestValid){bestValid=code.size();bestCode=code;save("best-valid.json",bestCode,0);log("larger_feasible_code");return true;}
        return false;
    }
    void verifyInternal(){vector<int>a(g.M);int actual=0;for(int w:code){assert(position[w]>=0);for(int c:g.d[w])if(a[c]++)++actual;}assert(actual==cost);assert(a==occ);}

    // Feasible destroy/repair representation. Counts refer to occupied received words.
    void fillGreedy(int heuristic,int sample,int banned=-1){
        heuristic%=4;
        vector<int> available;available.reserve(g.N);vector<int> where(g.N,-1),tight(g.N);
        auto remove=[&](int w){int k=where[w];if(k>=0){where[available.back()]=k;available[k]=available.back();available.pop_back();where[w]=-1;}};
        for(int w=0;w<g.N;++w)if(position[w]<0&&w!=banned){tight[w]=conflicts(w);if(!tight[w]){where[w]=available.size();available.push_back(w);}}
        int syndrome=random(g.n+1);
        while(!available.empty()&&timeLeft()){
            int choice=available[random(available.size())];double lowest=1e100;
            for(int k=0;k<sample;++k){
                int w=available[random(available.size())];double priority=unit()*.7;
                if(heuristic==0)priority+=g.d[w].size();
                if(heuristic==1){priority=unit();}
                if(heuristic==2){priority+=g.d[w].size()+3*(g.vt[w]!=syndrome);}
                if(heuristic==3){
                    int impact=0;
                    for(int c:g.d[w])for(int parent:g.parents[c])impact+=where[parent]>=0;
                    priority+=impact;
                }
                if(priority<lowest){lowest=priority;choice=w;}
            }
            remove(choice);add(choice);
            for(int c:g.d[choice])for(int w:g.parents[c])if(where[w]>=0)remove(w);
        }
        assert(cost==0);
    }

    void improveOneTwo(){
        // Standard (1,2) independent-set swaps: remove one word and add two.
        while(timeLeft()){
            vector<int> owner(g.M,-1);for(int w:code)for(int c:g.d[w])owner[c]=w;
            vector<vector<int>> groups(g.N);
            for(int w=0;w<g.N;++w)if(position[w]<0){
                int blocker=-1;bool multiple=false;
                for(int c:g.d[w])if(owner[c]>=0){if(blocker<0)blocker=owner[c];else if(blocker!=owner[c]){multiple=true;break;}}
                if(!multiple&&blocker>=0)groups[blocker].push_back(w);
            }
            vector<int> order=code;shuffle(order.begin(),order.end(),rng);bool changed=false;
            for(int old:order){
                auto &candidates=groups[old];shuffle(candidates.begin(),candidates.end(),rng);
                for(size_t i=0;i<candidates.size()&&!changed;++i)for(size_t j=i+1;j<candidates.size();++j){
                    int a=candidates[i],b=candidates[j];bool overlaps=false;size_t p=0,q=0;
                    while(p<g.d[a].size()&&q<g.d[b].size()){
                        if(g.d[a][p]==g.d[b][q]){overlaps=true;break;}
                        if(g.d[a][p]<g.d[b][q])++p;else ++q;
                    }
                    if(!overlaps){eraseIndex(position[old]);add(a);add(b);assert(cost==0);changed=true;break;}
                }
                if(changed)break;
            }
            if(!changed)break;
        }
    }

    void feasible(const string&policy,int heuristic,int sample,int maxShake){
        assert(cost==0);fillGreedy(heuristic,sample);if(heuristic>=4)improveOneTwo();record();vector<int> current=code;
        int failures=0;double lastLog=elapsed();
        while(timeLeft()){
            ++iteration;load(current);int baseline=code.size();
            int strength=policy=="greedy"?1:1+random(maxShake);
            if(policy=="adaptive")strength=min(maxShake,2+failures/30+random(4));
            int mode=policy=="restart"?3:random(3),removed=-1;
            if(policy=="restart"){load({});}
            else if(mode==0){
                for(int i=0;i<strength&&!code.empty();++i){int k=random(code.size());removed=code[k];eraseIndex(k);}
            }else if(mode==1){
                int center=code[random(code.size())],radius=1+random(5);
                for(int i=(int)code.size()-1;i>=0;--i)if(abs(g.weight[code[i]]-g.weight[center])<=radius&&unit()<.65){removed=code[i];eraseIndex(i);}
            }else{
                // Forced outsider / ejection: destroy all blockers, not just one weight layer.
                int outsider=random(g.N);
                if(position[outsider]<0){
                    vector<int> owners(g.M,-1);for(int w:code)for(int c:g.d[w])owners[c]=w;
                    vector<int> blockers;for(int c:g.d[outsider])if(owners[c]>=0)blockers.push_back(owners[c]);
                    sort(blockers.begin(),blockers.end());blockers.erase(unique(blockers.begin(),blockers.end()),blockers.end());
                    for(int w:blockers)eraseIndex(position[w]);add(outsider);
                }
            }
            fillGreedy(heuristic,sample,removed);if(heuristic>=4)improveOneTwo();record();
            int change=(int)code.size()-baseline;bool accept=change>0;
            if(policy!="greedy")accept=change>=0||unit()<exp(change/(policy=="adaptive"?1.8:.7));
            if(policy=="restart")accept=true;
            if(accept){current=code;if((int)code.size()>bestValid)failures=0;else ++failures;}else ++failures;
            if(policy=="adaptive"&&(int)current.size()<bestValid-8){current=bestCode;}
            if(policy=="adaptive"&&failures>300&&bestValid>0){current=bestCode;failures=0;}
            if(elapsed()-lastLog>2){verifyInternal();log("progress",",\"current_accepted_size\":"+to_string(current.size())+",\"failures\":"+to_string(failures));lastLog=elapsed();}
        }
    }

    void fixed(int target,const string&policy,int samples,int fullscanEvery){
        // Fixed cardinality permits invalid intermediate packings; cost=0 alone is valid.
        while((int)code.size()<target){
            int choice=-1,minimum=1000000;
            for(int j=0;j<g.N;++j){int w=random(g.N);if(position[w]>=0)continue;int c=conflicts(w);if(c<minimum){minimum=c;choice=w;}}
            assert(choice>=0);add(choice);
        }
        assert((int)code.size()==target);bestCost=1000000;record();
        int stagnant=0,penaltyUpdates=0;double lastLog=elapsed();
        while(timeLeft()&&bestCost>0){
            ++iteration;
            vector<int> conflicted;
            for(int i=0;i<(int)code.size();++i)if(responsibility(code[i]))conflicted.push_back(i);
            assert(!conflicted.empty());
            int chosen=conflicted[random(conflicted.size())];
            // Mixed focus can deliberately perturb the currently conflict-free backbone.
            if(policy!="conflict_only"&&unit()<.05)chosen=random(code.size());
            int old=code[chosen];eraseIndex(chosen);
            int newWord=-1,minimum=100000000;double tie=-1;bool bestAspiration=false;
            auto consider=[&](int w){
                if(w==old||position[w]>=0)return;
                int c=0,weighted=0;for(int child:g.d[w])if(occ[child]){++c;weighted+=penalties[child];}
                ++evaluations;
                if(tabu[w]>iteration&&cost+c>=bestCost)return;
                if(policy!="guided"&&policy!="breakout")weighted=c;
                bool aspiration=cost+c<bestCost;
                double t=unit();
                if((aspiration&&!bestAspiration)||(aspiration==bestAspiration&&(weighted<minimum||(weighted==minimum&&t>tie)))){minimum=weighted;newWord=w;tie=t;bestAspiration=aspiration;}
            };
            for(int bit=0;bit<g.n;++bit)consider(old^(1<<bit));
            for(int k=0;k<samples;++k){
                if(k%3==0)consider(random(g.N));
                else if(k%3==1)consider(old^(1<<random(g.n))^(1<<random(g.n)));
                else{
                    int i=random(g.n),j=random(g.n),b=random(2);
                    int shortened=(old&((1<<i)-1))|((old>>(i+1))<<i);
                    consider((shortened&((1<<j)-1))|((shortened>>j)<<(j+1))|(b<<j));
                }
            }
            if(fullscanEvery>0&&iteration%fullscanEvery==0)for(int w=0;w<g.N;++w)consider(w);
            if(newWord<0){newWord=old;minimum=conflicts(old);}
            add(newWord);tabu[old]=iteration+7+random(25);
            int previousBest=bestCost;record();
            if(bestCost<previousBest)stagnant=0;else ++stagnant;
            if((policy=="guided"||policy=="breakout")&&iteration%100==0){
                for(int c=0;c<g.M;++c)if(occ[c]>1)++penalties[c];
                ++penaltyUpdates;
            }
            if(policy=="guided"&&(cost>bestCost+25||stagnant>30000)){
                // Retain learned penalties while returning from an unproductive excursion.
                load(bestState);stagnant=0;
                log("bounded_excursion_restart");
            }
            // Explicit barrier crossing at the level of a collection of words.
            if(policy=="barrier"&&stagnant>2000){
                int edits=3+random(15);
                for(int j=0;j<edits;++j){int i=random(code.size());eraseIndex(i);int w=random(g.N);while(position[w]>=0)w=random(g.N);add(w);}
                log("barrier_crossing",",\"edits\":"+to_string(edits));stagnant=0;
            }
            if(elapsed()-lastLog>2){verifyInternal();log("progress",",\"stagnation\":"+to_string(stagnant)+",\"conflicted_words\":"+to_string(conflicted.size())+",\"penalty_updates\":"+to_string(penaltyUpdates));lastLog=elapsed();}
        }
    }

    void exchangeProgram(){
        // Outsiders are a feasible packing; parent words are restored automatically
        // whenever no selected outsider blocks them. Optimize gain |B|-|N_A(B)|.
        const auto parent=code;int baseSize=parent.size(),covered=0;
        vector<int> baseOwner(g.M,-1),cover(baseSize),owner(g.M,-1);
        for(int i=0;i<baseSize;++i)for(int c:g.d[parent[i]])baseOwner[c]=i;
        vector<vector<int>> labels(g.N),byLabel(baseSize);vector<int> outsiders,low;
        for(int w=0;w<g.N;++w)if(!initialMask[w]){
            auto &r=labels[w];for(int c:g.d[w])if(baseOwner[c]>=0)r.push_back(baseOwner[c]);
            sort(r.begin(),r.end());r.erase(unique(r.begin(),r.end()),r.end());
            for(int i:r)byLabel[i].push_back(w);outsiders.push_back(w);
            if(r.size()<=3)low.push_back(w);
        }
        load({});
        auto addB=[&](int w){add(w);assert(cost==0);for(int c:g.d[w])owner[c]=w;for(int i:labels[w])if(cover[i]++==0)++covered;};
        auto removeB=[&](int w){for(int c:g.d[w])owner[c]=-1;for(int i:labels[w])if(--cover[i]==0)--covered;eraseIndex(position[w]);};
        auto blockers=[&](int w){vector<int> b;for(int c:g.d[w])if(owner[c]>=0)b.push_back(owner[c]);sort(b.begin(),b.end());b.erase(unique(b.begin(),b.end()),b.end());return b;};
        auto fullCode=[&](){auto f=code;for(int i=0;i<baseSize;++i)if(!cover[i])f.push_back(parent[i]);return f;};
        set<pair<int,int>> exchangeArchive;double lastLog=elapsed();int stagnant=0,bestDistance=0;
        vector<int> delta(baseSize),touched,stamp(baseSize);int epoch=0;
        while(timeLeft()&&bestValid<baseSize+1){
            ++iteration;int choice=-1,chosenScore=-1000000,chosenSize=0;double priority=-1e100;vector<int> chosenBlockers;
            for(int sample=0;sample<384;++sample){
                int w;
                if(sample%3==0&&!low.empty())w=low[random(low.size())];
                else if(sample%3==1&&covered){int p=random(baseSize);while(!cover[p])p=random(baseSize);w=byLabel[p][random(byLabel[p].size())];}
                else w=outsiders[random(outsiders.size())];
                if(position[w]>=0||tabu[w]>iteration)continue;
                auto b=blockers(w);++epoch;touched.clear();
                auto change=[&](int p,int v){if(stamp[p]!=epoch){stamp[p]=epoch;delta[p]=0;touched.push_back(p);}delta[p]+=v;};
                for(int old:b)for(int p:labels[old])change(p,-1);
                for(int p:labels[w])change(p,1);
                int newCovered=covered;
                for(int p:touched)newCovered+=(cover[p]+delta[p]>0)-(cover[p]>0);
                int newSize=(int)code.size()-(int)b.size()+1,score=baseSize-newCovered+newSize;
                double u=max(1e-12,unit());double p=score+.035*newSize-.55*std::log(-std::log(u));
                if(score>bestValid)p+=1000; // Never trade a valid larger code for surrogate reward.
                if(p>priority){priority=p;choice=w;chosenScore=score;chosenSize=newSize;chosenBlockers=std::move(b);}
            }
            if(choice<0)continue;
            for(int w:chosenBlockers){removeB(w);tabu[w]=iteration+10+random(30);}addB(choice);
            assert((int)code.size()==chosenSize);assert(baseSize-covered+(int)code.size()==chosenScore);
            if(chosenScore>bestValid){bestValid=chosenScore;bestCode=fullCode();save("best-valid.json",bestCode,0);log("larger_exchange_code",",\"full_size\":"+to_string(chosenScore));stagnant=0;}
            else ++stagnant;
            if(chosenScore>=baseSize-16){
                auto bin=make_pair(baseSize-chosenScore,(int)code.size()/8);
                if(exchangeArchive.insert(bin).second)save("archive/exchange-loss"+to_string(bin.first)+"-d"+to_string(bin.second)+".json",fullCode(),0);
            }
            bestDistance=max(bestDistance,(int)code.size());
            if(stagnant>30000){
                // A whole outsider cluster is discarded; other restarts preserve far states in the archive.
                while(!code.empty())removeB(code.back());fill(tabu.begin(),tabu.end(),0);stagnant=0;
                log("exchange_restart");
            }
            if(elapsed()-lastLog>2){verifyInternal();log("exchange_progress",",\"full_size\":"+to_string(chosenScore)+",\"parents_displaced\":"+to_string(covered)+",\"best_distance\":"+to_string(bestDistance));lastLog=elapsed();}
        }
        log("exchange_finish",",\"full_size\":"+to_string(baseSize-covered+(int)code.size())+",\"parents_displaced\":"+to_string(covered)+",\"best_distance\":"+to_string(bestDistance));
    }
};

int main(int argc,char**argv){
    if(argc==4&&string(argv[1])=="audit"){
        Data g(stoi(argv[2]));ofstream f(argv[3],ios::binary);
        for(auto& row:g.d){uint32_t count=row.size();f.write((char*)&count,4);for(uint32_t v:row)f.write((char*)&v,4);}
        return 0;
    }
    // n seconds seed family policy input-ints output-dir target heuristic sample shake/fullscan
    if(argc<12){cerr<<"n seconds seed family policy seed-ints out target heuristic samples strength\n";return 2;}
    int n=stoi(argv[1]);assert(n>=6&&n<=16);Data g(n);
    Search s(g,stoi(argv[3]),stod(argv[2]),argv[7]);
    vector<int> input;ifstream f(argv[6]);int w;while(f>>w){assert(w>=0&&w<g.N);input.push_back(w);}
    s.load(input);s.initialMask.assign(g.N,0);for(int w:input)s.initialMask[w]=1;
    s.verifyInternal();s.record();s.log("start",",\"family\":\""+string(argv[4])+"\",\"policy\":\""+string(argv[5])+"\"");
    if(string(argv[4])=="exchange")s.exchangeProgram();
    else if(string(argv[4])=="feasible")s.feasible(argv[5],stoi(argv[9]),stoi(argv[10]),stoi(argv[11]));
    else s.fixed(stoi(argv[8]),argv[5],stoi(argv[10]),stoi(argv[11]));
    s.verifyInternal();s.log("finish");
    if(!s.bestCode.empty())s.save("best-valid.json",s.bestCode,0);
    if(!s.bestState.empty())s.save("best-state.json",s.bestState,s.bestCost);
}

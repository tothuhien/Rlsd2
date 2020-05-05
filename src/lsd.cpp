//    LSD - Least Square Dating for etimating substitution rate and divergence dates
//    Copyright (C) <2015> <Thu-Hien To>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "pr.h"
#include "options.h"
#include "readData.h"
#include "dating.h"
#include "utils.h"
#include "node.h"
#include "date.h"
#include "estimate_root.h"
#include "confidence_interval.h"
#include "outliers.h"
#include "lsd.h"
#include <ctime>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;
using namespace lsd;


int lsd2(Pr* &opt,vector<double>& rho,vector<double> &mrca)
{

    // initialise input and output streams
    InputOutputStream *io  = new InputOutputFile(opt);
    printInterface(*(io->outResult), opt);
    clock_t start = clock();
    double elapsed_time;
    if (io->inOutgroup){
        extrait_outgroup(io, opt);
    }
    ifstream gr(opt->rate.c_str());
    *(io->outTree1)<<"#NEXUS\n";
    *(io->outTree2)<<"#NEXUS\n";
    bool constraintConsistent=true;
    int s=0;
    double median_rate = opt->rho_min;
    int e;
    if (opt->partitionFile!="") {
        e = readPartitionFile(*(io->inPartition), opt);
        if (e==EXIT_FAILURE) return e;
        for (int i=0;i<=opt->ratePartition.size();i++){
            opt->givenRate.push_back(false);
        }
    }
    for (int y=1;y<=opt->nbData;y++){
        *(io->outResult)<<"\nTree "<<y<<" \n";
        cout<<"\nTREE "<<y<<endl;
        cout<<"*PROCESSING:"<<endl;
        cout<<"Reading the tree ... "<<endl;
        opt->init();
        Node** nodes;
        e = tree2data(*(io->inTree),opt,s,nodes);
        if (e==EXIT_FAILURE) return e;
        if (!opt->relative) {
            io->inDate->seekg(0);//BQM: set the stream position to 0 to read the same date file for each tree in the tree set
            e = readDateFile(*io->inDate, opt,nodes,constraintConsistent);
            if (e==EXIT_FAILURE) return e;
        }
        computeSuc_polytomy(opt,nodes);
        collapseUnInformativeBranches(opt,nodes);
        if (!opt->rooted){
            nodes = unrooted2rooted(opt,nodes);
        }
        if (opt->relative){
            for (int i=0;i<opt->nbINodes;i++) nodes[i]->removeConstraint();
            for (int i=opt->nbINodes;i<=opt->nbBranches;i++){
                nodes[i]->newPConstraint('p',opt->leaves);
            }
            Date* dateRoot = new Date("", 'p',opt->mrca,0,0);
            opt->internalConstraints.clear();
            opt->internalConstraints.push_back(dateRoot);
        }
        if (y==1){
            *(io->outTree1)<<"Begin trees;\n";
            *(io->outTree2)<<"Begin trees;\n";
        }
        if (opt->c == -1){
            double mbl;
            e = median_branch_lengths(opt,nodes,mbl);
            if (e==EXIT_FAILURE) return e;
            opt->b = max(mbl,10./opt->seqLength);
            if (opt->variance>0){
                cout<<"Parameter to adjust variances was set to "<<opt->b<<" (settable via option b)"<<endl;
                *(io->outResult)<<"Parameter to adjust variances was set to "<<opt->b<<" (settable via option b)\n";
            }
        } else {
            opt->b = opt->c;
        }
        computeVariance(opt,nodes);
        double br=0;
        if (opt->givenRate[0]){
            string line;
            if( getline(gr, line)) {
                vector<double> all_rates = read_double_from_line(line);
                opt->rho = all_rates[0];
                if (all_rates.size() > 1){
                    int i = 1;
                    while (i<=opt->ratePartition.size() && i<all_rates.size()){
                        opt->multiplierRate.push_back(all_rates[i]/opt->rho);
                        opt->givenRate[i] = true;
                        i++;
                    }
                }
            } else{
                std::ostringstream oss;
                oss<<"- There are less given rates than the number of given trees.\n";
                opt->warningMessage.push_back(oss.str());
                opt->givenRate[0] = false;
            }
        }
        constraintConsistent = initConstraint(opt, nodes);
        bool medianRateOK = true;
        if (opt->e>0) medianRateOK = calculateOutliers(opt,nodes,median_rate);
        else if (opt->minblen<0) medianRateOK = calculateMedianRate(opt,nodes,median_rate);
        imposeMinBlen(*(io->outResult),opt,nodes,median_rate,medianRateOK);
        if (!opt->constraint){//LD without constraints
            if (!constraintConsistent){
                ostringstream oss;
                oss<<"- There's conflict in the input temporal constraints.\n";
                opt->warningMessage.push_back(oss.str());
            }
            if (opt->estimate_root==""){//keep the given root
                cout<<"Dating without temporal constraints ..."<<endl;
                without_constraint_multirates(opt,nodes,true);
                output(br,y,opt,nodes,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                rho.push_back(opt->rho);
                mrca.push_back(nodes[0]->D);
            }
            else if (opt->estimate_root=="k"){
                cout<<"Estimating the root position on the branch defined by given root ..."<<endl;
                double br=0;
                vector<int>::iterator iter=nodes[0]->suc.begin();
                int s1=(*iter);
                iter++;
                int s2=(*iter);
                br=nodes[s1]->B+nodes[s2]->B;
                nodes[s1]->V=variance(opt,br);
                nodes[s2]->V=nodes[s1]->V;
                without_constraint_active_set_lambda_multirates(br,opt,nodes,true);
                output(br,y,opt,nodes,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                mrca.push_back(nodes[0]->D);
                rho.push_back(opt->rho);
            }
            else{//estimate the root
                int r;
                if (opt->estimate_root.compare("l")==0){//improve locally the root around the given root
                    cout<<"Estimating the root position locally around the given root ..."<<endl;
                    r=estimate_root_without_constraint_local_rooted(opt,nodes);
                }
                else{ //forget the given root and re-estimate the position of the root over all branhces
                    cout<<"Estimating the root position on all branches ..."<<endl;
                    r=estimate_root_without_constraint_rooted(opt,nodes);
                }
                if (r>0){
                    Node** nodes_new = cloneLeaves(opt,nodes,0);
                    vector<int>::iterator iter=nodes[0]->suc.begin();
                    int s1=(*iter);
                    iter++;
                    int s2=(*iter);
                    for (int i=opt->nbINodes; i<=opt->nbBranches; i++) {
                        nodes_new[i]->status=nodes[i]->status;
                    }
                    reroot_rootedtree(br,r,s1,s2,opt,nodes,nodes_new);
                    without_constraint_active_set_lambda_multirates(br,opt,nodes_new,true);
                    output(br,y,opt,nodes_new,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                    mrca.push_back(nodes_new[0]->D);
                    rho.push_back(opt->rho);
                    for (int i=0;i<opt->nbBranches+1;i++) delete nodes_new[i];
                    delete[] nodes_new;
                }
            }
        }
        else {//QPD with temporal constrains
            if (constraintConsistent || (opt->estimate_root!="" && opt->estimate_root!="k")){
                if (constraintConsistent){
                    if (opt->estimate_root==""){//keep the given root
                        if (constraintConsistent){
                            cout<<"Dating under temporal constraints mode ..."<<endl;
                            constraintConsistent = with_constraint_multirates(opt,nodes,true);
                        }
                        if (constraintConsistent) {
                            output(br,y,opt,nodes,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                            mrca.push_back(nodes[0]->D);
                            rho.push_back(opt->rho);
                        }
                        else{
                            cout<<"*WARNING: There's conflict in the input temporal constraints."<<endl;
                        }
                    }
                    else if (opt->estimate_root=="k"){
                        if (constraintConsistent){
                            cout<<"Estimating the root position on the branch defined by given outgroups ..."<<endl;
                            double br=0;
                            vector<int>::iterator iter=nodes[0]->suc.begin();
                            int s1=(*iter);
                            iter++;
                            int s2=(*iter);
                            br=nodes[s1]->B+nodes[s2]->B;
                            nodes[s1]->V=variance(opt,br);
                            nodes[s2]->V=nodes[s1]->V;
                            constraintConsistent = with_constraint_active_set_lambda_multirates(br,opt,nodes,true);
                        }
                        if (constraintConsistent) {
                            output(br,y,opt,nodes,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                            mrca.push_back(nodes[0]->D);
                            rho.push_back(opt->rho);
                        } else {
                            cout<<"*WARNING: There's conflict in the input temporal constraints."<<endl;
                        }
                    }
                    else{//estimate the root
                        int r;
                        if (opt->estimate_root.compare("l")==0){//improve locally the root around the given root
                            cout<<"Estimating the root position locally around the given root ..."<<endl;
                            r=estimate_root_with_constraint_local_rooted(opt,nodes);
                        }
                        else if (opt->estimate_root.compare("a")==0){
                            //forget the given root and re-estimate the position of the root over all branhces using fast method
                            cout<<"Estimating the root position on all branches using fast method ..."<<endl;
                            r=estimate_root_with_constraint_fast_rooted(opt,nodes);
                        }
                        else{ //forget the given root and re-estimate the position of the root over all branches
                            cout<<"Estimating the root position on all branches ..."<<endl;
                            r=estimate_root_with_constraint_rooted(opt,nodes);
                        }
                        Node** nodes_new = cloneLeaves(opt,nodes,0);
                        vector<int>::iterator iter=nodes[0]->suc.begin();
                        int s1=(*iter);
                        iter++;
                        int s2=(*iter);
                        for (int i=opt->nbINodes; i<=opt->nbBranches; i++) {
                            nodes_new[i]->status=nodes[i]->status;
                        }
                        reroot_rootedtree(br,r,s1,s2,opt,nodes,nodes_new);
                        with_constraint_active_set_lambda_multirates(br,opt,nodes_new,true);
                        output(br,y,opt,nodes_new,*(io->outResult),*(io->outTree1),*(io->outTree2),*(io->outTree3));
                        mrca.push_back(nodes_new[0]->D);
                        rho.push_back(opt->rho);
                        for (int i=0;i<opt->nbBranches+1;i++) delete nodes_new[i];
                        delete[] nodes_new;
                    }
                } else {
                    cout<<"*WARNING: There's conflict in the input temporal constraints."<<endl;
                }
            } else {
                cout<<"*WARNING: There's conflict in the input temporal constraints."<<endl;
            }
        }
        for (int i=0;i<=opt->nbBranches;i++) delete nodes[i];
        delete[] nodes;
    }
    *(io->outResult)<<"\n*********************************************************\n";
    *(io->outTree1)<<"End;\n";
    *(io->outTree2)<<"End;\n";
    gr.close();
    delete io;
    return EXIT_SUCCESS;
}

extern "C" SEXP Rlsd2(SEXP inputTree, SEXP inputDate, SEXP partitionFile, SEXP outFile, SEXP outGroup,
                     SEXP givenRate, SEXP seqLen, SEXP constraint , SEXP variance ,
                     SEXP confidenceIntervalSampling , SEXP rhoMin , SEXP estimate_root , SEXP b ,
                     SEXP q , SEXP mrca , SEXP leaves , SEXP verbose , SEXP keepOutgroup ,
                     SEXP nullblen , SEXP support ,  SEXP minblen ,SEXP  minblenL ,
                     SEXP roundTime, SEXP outDateFormat ,SEXP m ,SEXP e, SEXP nbData){
    int argc = 1;
    vector<char*> options;
    char* s;
    if (asLogical(constraint)){
        options.push_back((char*)"-c");
        argc++;
    }
    if (asLogical(verbose)){
        options.push_back((char*)"-j");
        argc++;
    }
    if (asLogical(keepOutgroup)){
        options.push_back((char*)"-k");
        argc++;
    }
    if (TYPEOF(inputDate) == STRSXP && STRING_PTR(inputDate)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(inputDate)[0]);
        options.push_back((char*)"-d");
        options.push_back(s);
        argc = argc + 2;
    }
    if (TYPEOF(partitionFile) == STRSXP && STRING_PTR(partitionFile)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(partitionFile)[0]);
        options.push_back((char*)"-p");
        options.push_back(s);
        argc = argc + 2;
    }
    if (TYPEOF(outFile) == STRSXP && STRING_PTR(outFile)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(outFile)[0]);
        options.push_back((char*)"-o");
        options.push_back(s);
        argc = argc + 2;
    }
    if (TYPEOF(outGroup) == STRSXP && STRING_PTR(outGroup)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(outGroup)[0]);
        options.push_back((char*)"-g");
        options.push_back(s);
        argc = argc + 2;
    }
    if (TYPEOF(givenRate) == STRSXP && STRING_PTR(givenRate)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(givenRate)[0]);
        options.push_back((char*)"-w");
        options.push_back(s);
        argc = argc + 2;
    }
    if (TYPEOF(estimate_root) == STRSXP && STRING_PTR(estimate_root)[0]!=NA_STRING){
        s = (char*)CHAR(STRING_PTR(estimate_root)[0]);
        options.push_back((char*)"-r");
        options.push_back(s);
        argc = argc + 2;
    }
    double Q = asReal(q);
    if (asReal(q)>=0){
        options.push_back((char*)"-q");
        s = new char[size_t(Q)];
        sprintf(s, "%e", Q);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(estimate_root) == REALSXP && !ISNA(REAL(nullblen)[0])){
        double NULLBLEN = asReal(nullblen);
        options.push_back((char*)"-l");
        s = new char[size_t(NULLBLEN)];
        sprintf(s, "%e", NULLBLEN);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(support) == REALSXP && !ISNA(REAL(support)[0])){
        double SUPPORT = asReal(support);
        options.push_back((char*)"-S");
        s = new char[size_t(SUPPORT)];
        sprintf(s, "%e", SUPPORT);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(minblen) == REALSXP && !ISNA(REAL(minblen)[0])){
        double MINBLEN = asReal(minblen);
        options.push_back((char*)"-u");
        s = new char[size_t(MINBLEN)];
        sprintf(s, "%e", MINBLEN);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(minblenL) == REALSXP && !ISNA(REAL(minblenL)[0])){
        double MINBLENL = asReal(minblenL);
        options.push_back((char*)"-U");
        s = new char[size_t(MINBLENL)];
        sprintf(s, "%e", MINBLENL);
        options.push_back(s);
        argc = argc+2;
    }
    if (INTEGER(confidenceIntervalSampling)[0]!=NA_INTEGER){
        int NBS = asInteger(confidenceIntervalSampling);
        options.push_back((char*)"-f");
        s = new char[size_t(NBS)];
        sprintf(s, "%d", NBS);
        options.push_back(s);
        argc = argc+2;
    }
    if (INTEGER(outDateFormat)[0]!=NA_INTEGER){
        int D = asInteger(outDateFormat);
        options.push_back((char*)"-D");
        s = new char[size_t(D)];
        sprintf(s, "%d", D);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(b) == REALSXP && !ISNA(REAL(b)[0])){
        double B = asReal(b);
        options.push_back((char*)"-b");
        s = new char[size_t(B)];
        sprintf(s, "%e", B);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(mrca) == REALSXP && !ISNA(REAL(mrca)[0])){
        double A = asReal(mrca);
        options.push_back((char*)"-a");
        s = new char[size_t(A)];
        sprintf(s, "%e", A);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(leaves) == REALSXP && !ISNA(REAL(leaves)[0])){
        double Z = asReal(leaves);
        options.push_back((char*)"-z");
        s = new char[size_t(Z)];
        sprintf(s, "%e", Z);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(roundTime) == REALSXP && !ISNA(REAL(roundTime)[0])){
        double ROUND_TIME = asReal(roundTime);
        options.push_back((char*)"-R");
        s = new char[size_t(ROUND_TIME)];
        sprintf(s, "%e", ROUND_TIME);
        options.push_back(s);
        argc = argc+2;
    }
    if (TYPEOF(e) == REALSXP && !ISNA(REAL(e)[0])){
        double E = asReal(e);
        options.push_back((char*)"-e");
        s = new char[size_t(E)];
        sprintf(s, "%e", E);
        options.push_back(s);
        argc = argc+2;
    }
    int var = asInteger(variance);
    if (var != 0){
        options.push_back((char*)"-v");
        s = new char[size_t(var)];
        sprintf(s, "%d", var);
        options.push_back(s);
        argc = argc+2;
    }
    argc = argc + 10;
    char** argv = new char*[argc];
    argv[1] = (char*)"-i";
    argv[2] = (char*)CHAR(STRING_PTR(inputTree)[0]);
    argv[3] = (char*)"-s";
    int sl = asInteger(seqLen);
    argv[4] = new char[size_t(sl)];
    sprintf(argv[4], "%d", sl);
    argv[5] = (char*)"-t";
    argv[6] = new char[size_t(asReal(rhoMin))];
    sprintf(argv[6], "%e", asReal(rhoMin));
    argv[7] = (char*)"-m";
    argv[8] = new char[size_t(asInteger(m))];
    sprintf(argv[8], "%d", asInteger(m));
    argv[9] = (char*)"-n";
    argv[10] = new char[size_t(asInteger(nbData))];
    sprintf(argv[10], "%d", asInteger(nbData));
    for (int i = 11; i<argc;i++){
        argv[i] = options[i-11];
    }
    /*for (int i=0;i<options.size();i++){
      cout<<options[i]<<endl;
    }*/
    Pr* opt = getOptions( argc, argv);
    if (opt==NULL){
        return R_NilValue;
    }
    vector<double> rates;
    vector<double> rootDate;
    int r = lsd2(opt,rates,rootDate);
    if (r==EXIT_SUCCESS){
        SEXP rho = PROTECT(allocVector(REALSXP, rates.size()));
        SEXP root = PROTECT(allocVector(REALSXP, rootDate.size()));
        for (int i=0;i<rates.size();i++){
            REAL(rho)[i] = rates[i];
            REAL(root)[i] = rootDate[i];
        }
        SEXP res = PROTECT(allocVector(VECSXP, 2));
        SET_VECTOR_ELT(res, 0, rho);
        SET_VECTOR_ELT(res, 1, root);
        UNPROTECT(3);
        delete opt;
        return res;
    } else {
        return R_NilValue;
    }
}
//R CMD SHLIB confidence_interval.cpp dating.cpp estimate_root.cpp options.cpp outliers.cpp readData.cpp utils.cpp lsd.cpp -o Rlsd2.so
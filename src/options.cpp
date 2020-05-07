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
#include "options.h"

Pr* getOptions( int argc, char** argv )
{
    return getCommandLine( argc, argv );
}

Pr* getCommandLine( int argc, char** argv)
{
    const string VERSION="v1.6.6";
    Pr* opt = new Pr();
    int c;
    string s;
    bool iflag = false,
    dflag = false,
    flagA=false,
    flagZ=false,
    sflag=false,
    fflag=false,
    lflag=false,
    uflag=false,
    vflag=false,
    validDate = true;
    optind = 1;
    while ( (c = getopt(argc, argv, ":i:d:D:o:s:n:g:r:v:ct:w:b:a:z:f:kje:m:p:q:u:l:U:R:S:")) != -1 )
    {
        switch (c)
        {
            case 'i':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 ){
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                    return NULL;
                }
#endif
                opt->inFile = optarg;
                iflag = true;
                break;
            case 'd':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 ){
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                    return NULL;
                }
#endif
                opt->inDateFile = optarg;
                dflag = true;
                break;
            case 'D':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -D must be an integer.\n");
                    return NULL;
                }
                opt->outDateFormat = atoi(optarg);
                if (opt->outDateFormat !=1 && opt->outDateFormat != 2){
                    myExit("Argument of option -D must be either 1 (date as real) or 2 (date as YY-MM-DD).\n");
                    return NULL;
                }
                break;
            case 'p':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 ){
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                    return NULL;
                }
#endif
                opt->partitionFile = optarg;
                break;
            case 'o':
                opt->outFile = optarg;
                break;
            case 's':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -s must be an integer.\n");
                    return NULL;
                }
                opt->seqLength = atoi(optarg);
                if( opt->seqLength<1 ){
                    myExit("Argument of option -s must be strictly positive.\n");
                    return NULL;
                }
                sflag = true;
                break;
            case 'n':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -n must be an integer.\n");
                    return NULL;
                }
                opt->nbData = atoi(optarg);
                if( opt->nbData<1 ){
                    myExit("Argument of option -n must be strictly positive.\n");
                    return NULL;
                }
                break;
            case 'g':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 ){
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                    return NULL;
                }
#endif
                opt->fnOutgroup = optarg;
                break;
            case 'k':
                if (opt->fnOutgroup!=""){
                	opt->keepOutgroup=true;
                	opt->estimate_root="k";
                }
                break;
            case 'r':
                opt->estimate_root = optarg;
                if (opt->estimate_root.compare("l")!=0 && opt->estimate_root.compare("a")!=0 && opt->estimate_root.compare("as")!=0 && opt->estimate_root.compare("k")!=0){
                    myExit("Argument of option -r must be either \"l\" or \"a\" or \"as\".\n");
                    return NULL;
                }
                break;
            case 'v':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -v must be either 1 or 2.\n");
                    return NULL;
                }
                opt->variance = atoi(optarg);
                if (opt->variance!=1 && opt->variance!=2){
                    myExit("Argument of option -v must be either 1 (for using orginal branches to compute variance) or 2 (LSD will be run twice, the second time uses the variances based on the estimated branch lengths of the first time).\n");
                    return NULL;
                }
                vflag = true;
                break;
            case 'c':
                opt->constraint = true;
                break;
            case 'b':
                if( !isReal(optarg) ){
                    myExit("Argument of option -b must be a real.\n");
                    return NULL;
                }
                opt->c = atof( optarg );
                if (opt->c<=0 || opt->c>1){
                    myExit("Argument of option -b must be a positive number samller than 1, see the help page of lsd2 -h for more information.\n");
                    return NULL;
                }
                break;
            case 'e':
                if( !isReal(optarg) ){
                    myExit("Argument of option -e must be a real.\n");
                    return NULL;
                }
                opt->e = atof(optarg);
                if (opt->e<0){
                    std::ostringstream oss;
                    oss<<"- The specified argument of option e was negative, so outliers detection option was not processed\n";
                    opt->warningMessage.push_back(oss.str());
                }
                break;
            case 'm':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -m must be an integer.\n");
                    return NULL;
                }
                opt->m = atof( optarg );
                if (opt->m<2) {
                    myExit("Argument of option -m must be >= 2.\n");
                    return NULL;
                }
                break;
            case 't':
                if( !isReal(optarg) ){
                    myExit("Argument of option -t must be a real.\n");
                    return NULL;
                }
                opt->rho_min = atof(optarg);
                if (opt->rho_min<=0){
                    myExit("Argument of option -t must be strictly positive.\n");
                    return NULL;
                }
                break;
            case 'q':
                if (!isReal(optarg)){
                    myExit("Argument of option -q must be a real.\n");
                    return NULL;
                }
                opt->q = atof(optarg);
                if (opt->q<0){
                    myExit("Argument of option -q could not be negative.\n");
                    return NULL;
                }
                break;
            case 'w':
                if( access( optarg, R_OK )!=0 ){
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                    return NULL;
                }
                opt->rate = optarg;
                opt->givenRate[0] = true;
                break;
            case 'a':
                if( isReal(optarg) ){
                    opt->mrca=atof(optarg);
                    validDate=true;
                }
                else{
                    validDate=readDateFromString(optarg,opt->mrca);
                    opt->inDateFormat=2;
                }
                if (!validDate){
                    myExit("Argument of option -a must be a real or a date format year-month-day.\n");
                    return NULL;
                }
                flagA=true;
                break;
            case 'z':
                if( isReal(optarg) ){
                    opt->leaves=atof(optarg);
                    validDate=true;
                }
                else{
                    validDate=readDateFromString(optarg,opt->leaves);
                    opt->inDateFormat=2;
                }
                if (!validDate){
                    myExit("Argument of option -z must be a real or a date format year-month-day.\n");
                    return NULL;
                }
                flagZ=true;
                break;
            case 'j':
                opt->verbose = true;
                break;
            case 'f':
                if( !isInteger(optarg) ){
                    myExit("Argument of option -f must be the number of samplings to compute confidence intervals, e.g. 100 ...\n");
                    return NULL;
                }
                opt->nbSampling = atoi( optarg );
                if (opt->nbSampling<0){
                    myExit("Argument of option -f must be a positive integer.\n");
                    return NULL;
                }
                fflag = true;
                opt->ci=true;
                break;
            case 'u':
                if ( !isReal(optarg)){
                    myExit("Argument of option -u must be a real\n");
                    return NULL;
                }
                opt->minblen = atof(optarg);
                if (opt->minblen<0){
                    myExit("Argument of option -u must be >= 0\n");
                    return NULL;
                }
                uflag = true;
                break;
            case 'U':
                if ( !isReal(optarg)){
                    myExit("Argument of option -U must be a real\n");
                    return NULL;
                }
                opt->minblenL = atof(optarg);
                if (opt->minblenL<0){
                    myExit("Argument of option -U must be >= 0\n");
                    return NULL;
                }
                uflag = true;
                break;
            case 'l':
                if ( !isReal(optarg)){
                    myExit("Argument of option -l must be a real\n");
                    return NULL;
                }
                opt->nullblen = atof(optarg);
                lflag = true;
                break;
            case 'R':
                if ( !isReal(optarg)){
                    myExit("Argument of option -R must be a real\n");
                    return NULL;
                }
                opt->round_time = atof(optarg);
                if (opt->round_time<=0){
                    myExit("Argument of option -R must be positive\n");
                    return NULL;
                }
                break;
            case 'S':
                if ( !isReal(optarg)){
                    myExit("Argument of option -S must be a real\n");
                    return NULL;
                }
                opt->support = atof(optarg);
                if (opt->round_time<0){
                    myExit("Argument of option -S must be positive\n");
                    return NULL;
                }
                break;
            case '?':
                myExit("Unrecognized option: -%c\n", optopt);
                return NULL;
            case ':':
                if (optopt=='v'){
                    myExit("Argument of option -v must be either 1 (for using orginal branches to compute variance) or 2 (for using branches estimated by LSD to compute variance, i.e LSD will be run 2 times: the first time is just used to compute the variance).\n");
                    return NULL;
                }
                else {
                    myExit("Option -%c requires an operand\n", optopt );
                    return NULL;
                }
            default:
                myExit("?? getopt returned character code 0%o ??\n", c );
                return NULL;
        }
    }
    if( !(iflag) ){
        myExit("Argument -i is necessary to continue.\n");
        return NULL;
    }
    if (!dflag && (flagA && !flagZ)){
        myExit("The input date file is not provided, so option -z is required together with option -a to estimate relative dates.\n");
        return NULL;
    }
    if (!dflag && (!flagA && flagZ)){
        myExit("The input date file is not provided, so option -a is required together with option -z to estimate relative dates.\n");
        return NULL;
    }
    if ( !(sflag) && (fflag || !uflag || !lflag || vflag)){
        myExit("Argument -s is necessary to continue.\n");
        return NULL;
    }
    if (!dflag && (!flagA && !flagZ)){
        opt->relative=true;
        opt->mrca=0;
        opt->leaves=1;
    }
    if (!dflag && (flagA && flagZ)){
        if (opt->inDateFormat==2 && opt->outDateFormat==0) opt->outDateFormat=2;
        if (opt->mrca >= opt->leaves){
            myExit("The root date must be strictly smaller than the tips date.\n");
            return NULL;
        }
        opt->relative=true;
    }

    if (dflag)
        opt->relative=false;

    if (!opt->constraint && opt->estimate_root.compare("as")==0){
        cout<<"The non constrained mode is chosen, so the \"as\" method for rooting function is the same as the \"a\" method."<<endl;
        opt->estimate_root="a";
    }
    if( opt->outFile=="") opt->outFile = opt->inFile + ".result";
    opt->treeFile1=opt->outFile+".nexus";
    opt->treeFile2=opt->outFile+".date.nexus";
    opt->treeFile3=opt->outFile+".nwk";
    if (std::isnan(opt->nullblen)){
        opt->nullblen = 0.5/opt->seqLength;
    }

    return opt;
}



void printInterface(ostream& in, Pr* opt)
{
    const string VERSION = "v1.6.6";

    in<<"\nLEAST-SQUARE METHODS TO ESTIMATE RATES AND DATES - "<<VERSION<<" \n\n";
    in<<"\nInput files:\n";
    in<<"  i                                               Input tree file : "<<opt->inFile.c_str()<<"\n";
    if (opt->relative==true){
        ostringstream tMRCA,tLeaves;
        if (opt->inDateFormat==2){
            tMRCA<<realToYearMonthDay(opt->mrca);
            tLeaves<<realToYearMonthDay(opt->leaves);
        } else {
            tMRCA<<opt->mrca;
            tLeaves<<opt->leaves;
        }
        in<<"  d                                       Estimate relative dates : mrca date = "<<tMRCA.str()<<", tips date = "<<tLeaves.str()<<"\n";
    }
    else
        in<<"  d                                               Input date file : "<<opt->inDateFile.c_str()<<"\n";
    in<<"  p                                                Partition file : ";
    if (opt->partitionFile=="")        in<<"No\n";
    else in<<opt->partitionFile.c_str()<<"\n";
    in<<"Output file:\n";
    in<<"  o                                                  Output file  : "<<opt->outFile.c_str()<<"\n";
    in<<"Parameters:\n";
    in<<"  c                                              With constraints : ";
    if (!opt->constraint) in<<"No\n";
    else {
        in<<"Yes\n";
    }
    in<<"  t                                      Lower bound for the rate : "<<opt->rho_min<<"\n";
    in<<"  v                                                With variances : ";
    if (opt->variance==0) in<<"No\n";
    else {
        if (opt->variance==1) in<<"Yes, use variances based on input branch lengths\n";
        else if (opt->variance==2) in<<"Yes, use variances based on estimated branch lengths\n";
        in<<"  b                              Adjusted parameter for variances : ";
        if (opt->c==-1){
            in<<"To estimate\n";
        } else{
            in<<opt->c<<"\n";
        }
    }
    in<<"  r                                             Estimate the root : ";
    if (opt->estimate_root==""){
        in<<"No\n";
    }
    else if (opt->estimate_root.compare("l")==0){
        in<<"Around the given root\n";
    }
    else if (opt->estimate_root.compare("a")==0 && opt->constraint){
        in<<"Use fast method to search on all branches\n";
    }
    else if (opt->estimate_root.compare("a")==0 && !opt->constraint){
        in<<"Search on all branches\n";
    }
    else if (opt->estimate_root.compare("as")==0){
        in<<"Use constrained mode on all branches\n";
    }
    in<<"  w                                       Given substitution rate : ";
    if (opt->rate=="") in<<"No\n";
    else in<<opt->rate.c_str()<<"\n";
    if (opt->fnOutgroup=="")
        in<<"  g                                               Given outgroups : No\n";
    else {
        in<<"  g                                       File contains outgroups : "<<opt->fnOutgroup.c_str()<<"\n";
        if (opt->keepOutgroup) {
            in<<"  k                         Keep outgroups in the estimating tree :  Yes\n";
        }
        else{
            in<<"  k                         Keep outgroups in the estimating tree : No\n";
        }
    }
    in<<"  n                                             Multiple data set : ";
    if( opt->nbData< 2 )
        in<<"No\n";
    else
        in<<"Yes, "<<opt->nbData<<" data sets\n";
    in<<"  f                                  Compute confidence intervals : ";
    if (opt->ci){
        in<<"Yes, sampling "<<opt->nbSampling<<" times\n";
        in<<"  q                  Standard deviation of lognormal relaxed clock: "<<opt->q<<" (for computing confidence intervals)\n";
    }
    else
        in<<"No\n";
    if (opt->nullblen==-1 || opt->ci){
        in<<"  s                                               Sequence length : "<<opt->seqLength<<"\n";
    }
    in<<"  e                                          Exclude outlier tips : ";
    if (opt->e>0){
        in<<"Yes, detect and exclude outliers from the analysis\n";
        in<<"  m                   Number of sampling nodes to detect outliers : "<<opt->m<<"\n";
        in<<"  e                       The Zscore threshold to detect outliers : "<<opt->e<<"\n";
    }
    else
        in<<"No\n";
    in<<"  u                Min internal branch length of time scaled tree : ";
    if (opt->minblen==-1){
        in<<"To estimate\n";
        in<<"  R     Rounding number for min branch length of time scaled tree : ";
        if (opt->round_time>0) in<<opt->round_time<<"\n";
        else in<<"To guess\n";
    } else {
        in<<opt->minblen<<"\n";
    }
    in<<"  U                Min external branch length of time scaled tree : ";
    if (opt->minblenL==-1){
        if (opt->minblen>0) {
            opt->minblenL = opt->minblen;
            in<<opt->minblenL<<"\n";
        } else {
            in<<"To estimate\n";
        }
    } else {
        in<<opt->minblenL<<"\n";
    }
    in<<"  l                        Collapsed internal branch length limit : ";
    if (opt->nullblen==-1){
        in<<0.5/opt->seqLength<<"\n";
    } else {
        in<<opt->nullblen<<"\n";
    }
    in<<"  D                                            Output date format :";
    if (opt->outDateFormat==0){
        in<<" Based on input date format\n";
    }
    if (opt->outDateFormat==1){
        in<<" Real number\n";
    }
    if (opt->outDateFormat==2){
        in<<" Year-Month-Day\n";
    }
    in<<"\n  y to accept or type a letter to change an option (x = Exit) ";
}




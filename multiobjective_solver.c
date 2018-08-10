/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#ifdef CPLEX
    #include "cplex.h"
#else
    #include <glpk.h>
#endif

#include "point_class.h"
#include "sydneys_class.h"
#include "simplex_class.h"
#include "problem_class.h"
#include "multiobjective_solver.h"
#include "SimplexStore.h"

bool DEBUG = false;
bool SCAN_FOR_REPEATS = false;
bool SAVE_POINTS = true;
bool SCAN_FOR_NEGATIVE_NORMAL = false;
double EPSILON = .000000001;
bool DIM_REDUCE = true;
bool SCALE_OBJECTIVES = false;

int main(int argc, char **argv)
{
	/* initialize Cplex environment *********************************/
	int status = 0;
	string phrase = "";
	string phrase2 = "";
	string phraseCopy = "";
	string filename1 = "", filename2 = "";
	string temp = "";
	MultiobjectiveProblem myProb;
	bool minProb = true;
	#ifdef CPLEX
	    CPXENVptr  env = NULL;
        CPXLPptr   lp = NULL;
        CPXLPptr   lp2 = NULL;
    #else 
        glp_prob *lp;
        glp_prob *lp2;
    #endif
    ifstream fin;
    ofstream fout;
    bool done = false;
	
	#ifdef CPLEX
      	env = CPXopenCPLEX(&status);

      	if(env==NULL)
      	{
        		char errmsg[1024];
        		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
        		CPXgeterrorstring (env, status, errmsg);
        		printf ("%s", errmsg);
        		exit(0);
      	}
  	#endif
  	
	/******************************************************************/

  	/************* Set to 1 thread **********************************/
  	
  	#ifdef CPLEX
	    status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	    if ( status ) {
		    printf ("Failure to set threads to 1, error %d.\n",status);
		    exit(0);
	    }
	    
/*	    status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);*/
/*	    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);*/
/*	    if ( status ) {*/
/*		    printf ("Failure to set threads to 1, error %d.\n",status);*/
/*		    exit(0);*/
/*	    }*/
	#endif
  	
  	/******************************************************************/
  	
  	/************* Get Number of Objectives ***************************/
  	
	myProb.SetNumObj(atoi(argv[1]));
  	
  	/******************************************************************/
  	
  	
  	/************* Set any desired CPLEX parameters here **************/
  	
  	#ifdef CPLEX
	    status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_AbsGap, 0.0);
	    if ( status ) {
		    printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		    exit(0);
	    }
	    status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_RelGap, 0.0);
	    if ( status ) {
		    printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		    exit(0);
	    }
	
	    myProb.SetEnv(env);
	#else
	    glp_term_out(GLP_OFF); // terminal output is on, but I will want it off (on is default)
	#endif
	
//	status = CPXsetincumbentcallbackfunc(env, userincumbent, NULL);
  	
  	/******************************************************************/
  	
  	/** Create test Problems and read them in from the model files ****/
    
    phrase = argv[2];
    phrase = phrase.substr(phrase.find("."));
    filename1 = argv[2];
    filename2 = "temp0.mop";
    
    if(phrase[1] == 'm') // assume we are working with a single *.mop file
    {
        #ifndef CPLEX
/*        cout << "modifying file to not contain OBJSENSE" << endl;*/
            fin.open(filename1);
	        fout.open(filename2);
            while (getline(fin, phrase))
            {
                phrase2 = regex_replace(phrase, regex("^ +"), "");
                if((phrase2[0] == 'O' || phrase2[0] == 'o') && (phrase2[1] == 'B' || phrase2[1] == 'b') && (phrase2[2] == 'J' || phrase2[2] == 'j'))
                {
                    getline(fin, phrase);
                    phrase2 = regex_replace(phrase, regex("^ +"), "");
                    if((phrase2[0] == 'M' || phrase2[0] == 'm') && (phrase2[1] == 'A' || phrase2[1] == 'a') && (phrase2[2] == 'X' || phrase2[2] == 'x')) minProb = false;
                }
                else
                {
                    fout << phrase << endl;
                }
            }
            fin.close();
            fout.close();
            filename1 = filename2;
            filename2[4]++;
        #endif
        
        for(int i = 0; i < myProb.GetNumObj(); i++)
	    {
	        if(i != 0) 
	        {
	            done = false;
	            fin.open(filename1);
	            fout.open(filename2);
	            while (getline(fin, phrase))
                {
                    phrase2 = regex_replace(phrase, regex("^ +"), "");
/*                    cout << phrase2 << endl;*/
                    if(!done && ((phrase2[0] == 'N' || phrase2[0] == 'n') && (phrase2[1] == ' ' || phrase2[1] == '\t'))) 
                    {
/*                        cout << phrase2 << endl;*/
                        temp = phrase;
                        done = true;
                    }
                    else
                    {
                        if(temp.length() && (phrase2[0] != 'N' && phrase2[0] != 'n' && phrase2[0] != 'L' && phrase2[0] != 'l' && phrase2[0] != 'G' && phrase2[0] != 'g' && phrase2[0] != 'E' && phrase2[0] != 'e'))
                        {
/*                            cout << "here" << endl;*/
                            fout << temp << endl;
                            temp = "";
                        }
/*                        cout << phrase << endl;*/
                        fout << phrase << endl;
                    }
                }
                fin.close();
                fout.close();
                filename1 = filename2;
                filename2[4]++;
	        }
	        #ifdef CPLEX
	          	lp = CPXcreateprob (env,&status,filename1.c_str());
	          	if(lp==NULL) 
	          	{
	            		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", i+1, status);
	            		exit(0);
	            }
	            	
	            status = CPXreadcopyprob(env,lp,filename1.c_str(),NULL);
	          	if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(0);
	            }
	        #else
	            lp = glp_create_prob();
	            status = glp_read_mps(lp, GLP_MPS_FILE, NULL, filename1.c_str());
	            if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(0);
	            }
	            if(minProb) glp_set_obj_dir(lp, GLP_MIN);
	            else glp_set_obj_dir(lp, GLP_MAX);
	        #endif
	        
	        myProb.AddLP(lp);

        }
        filename2[4] = '0';
        
        #ifdef CPLEX
            for(int i = 0; i < myProb.GetNumObj() - 1; i++)
            {
                remove(filename2.c_str());
                filename2[4]++;
            }
        #else
            for(int i = 0; i < myProb.GetNumObj(); i++)
            {
                remove(filename2.c_str());
                filename2[4]++;
            }
        #endif
    }
    else if(phrase[1] == 'v') //assume we are working with a *.vlp file
    {
        vector< vector<double> > objCoefs;
        vector< vector<double> > conCoefs;
        vector<char> rowTypes;
        vector<double> rowLB;
        vector<double> rowUB;
        vector<char> colTypes;
        vector<double> colLB;
        vector<double> colUB;
        vector<int> indices;
        int rows = 0;
        int cols = 0;
        int objs = 0;
        int intVal = 0;
        int intVal2 = 0;
        int intZ = 0;
        double negInf = -100000000000000000000.;
        double posInf = 100000000000000000000.;
        
        cout << "VLP file type detected." << endl;
        
        fin.open(filename1);
        while(getline(fin, phrase)) //read the problem data from the *.vlp file
        {
            phrase2 = regex_replace(phrase, regex("^ +"), "");
            phrase2 = regex_replace(phrase2, regex("^\t+"), "");
/*            cout << phrase2 << endl;*/
            if(phrase2[0] == 'p' || phrase2[0] == 'P')
            {
                phrase2 = phrase2.substr(phrase2.find(" "));
                phrase2 = regex_replace(phrase2, regex("^ +"), "");
                phrase2 = regex_replace(phrase2, regex("^\t+"), "");
/*                cout << phrase2 << endl;*/
                phrase2 = phrase2.substr(phrase2.find(" "));
                phrase2 = regex_replace(phrase2, regex("^ +"), "");
                phrase2 = regex_replace(phrase2, regex("^\t+"), "");
/*                cout << phrase2 << endl;*/
                if((phrase2[0] == 'M' || phrase2[0] == 'm') && (phrase2[1] == 'A' || phrase2[1] == 'a') && (phrase2[2] == 'X' || phrase2[2] == 'x')) minProb = false;
                phrase2 = phrase2.substr(phrase2.find(" "));
                phrase2 = regex_replace(phrase2, regex("^ +"), "");
                phrase2 = regex_replace(phrase2, regex("^\t+"), "");
/*                cout << phrase2 << endl;*/
                stringstream ss(phrase2);
                ss >> rows >> cols >> intVal >> objs;
                cout << "Number of rows: " << rows << endl;
                cout << "Number of cols: " << cols << endl;
                cout << "Number of objs: " << objs << endl;
                conCoefs.resize(rows);
                objCoefs.resize(objs);
                for(int i = 0; i < rows; i++)
                {
                    #ifdef CPLEX
                        conCoefs[i].resize(cols,0.);
                    #else
                        conCoefs[i].resize(cols+1,0.);
                    #endif
                }
                for(int i = 0; i < objs; i++)
                {
                    #ifdef CPLEX
                        objCoefs[i].resize(cols,0.);
                    #else
                        objCoefs[i].resize(cols+1,0.);
                    #endif
                }
                rowLB.resize(rows + 1,negInf);
                rowUB.resize(rows + 1,posInf);
                rowTypes.resize(rows + 1);
                colLB.resize(cols + 1,negInf);
                colUB.resize(cols + 1,posInf);
                colTypes.resize(cols + 1);
            }
            else if(phrase2[0] == 'i' || phrase2[0] == 'I')
            {
                phrase2 = phrase2.substr(phrase2.find(" "));
                stringstream ss(phrase2);
                ss >> intVal;
                ss >> rowTypes[intVal];
/*                cout << "Setting row " << intVal << " to " << rowTypes[intVal] << " type." << endl;*/
                if(rowTypes[intVal] == 'l') ss >> rowLB[intVal];
                else if(rowTypes[intVal] == 'u') ss >> rowUB[intVal];
                else if(rowTypes[intVal] == 'd') ss >> rowLB[intVal] >> rowUB[intVal];
                else if(rowTypes[intVal] == 's')
                {
                    ss >> rowLB[intVal];
                    rowUB[intVal] = rowLB[intVal];
                }
            }
            else if(phrase2[0] == 'j' || phrase2[0] == 'J')
            {
                phrase2 = phrase2.substr(phrase2.find(" "));
                stringstream ss(phrase2);
                ss >> intVal;
                ss >> colTypes[intVal];
/*                cout << "Setting col " << intVal << " to " << colTypes[intVal] << " type." << endl;*/
                if(colTypes[intVal] == 'l') ss >> colLB[intVal];
                else if(colTypes[intVal] == 'u') ss >> colUB[intVal];
                else if(colTypes[intVal] == 'd') ss >> colLB[intVal] >> colUB[intVal];
                else if(colTypes[intVal] == 's')
                {
                    ss >> colLB[intVal];
                    colUB[intVal] = colLB[intVal];
                }
            }
            else if(phrase2[0] == 'a' || phrase2[0] == 'A')
            {
                phrase2 = phrase2.substr(phrase2.find(" "));
                stringstream ss(phrase2);
                ss >> intVal >> intVal2;
                #ifdef CPLEX
                    ss >> conCoefs[intVal - 1][intVal2 - 1];
                #else
                    ss >> conCoefs[intVal - 1][intVal2];
                #endif
            }
            else if(phrase2[0] == 'o' || phrase2[0] == 'O')
            {
                phrase2 = phrase2.substr(phrase2.find(" "));
                stringstream ss(phrase2);
                ss >> intVal >> intVal2;
                #ifdef CPLEX
                    ss >> objCoefs[intVal - 1][intVal2 - 1];
                #else
                    ss >> objCoefs[intVal - 1][intVal2];
                #endif
            }
            else if(phrase2[0] == 'k' || phrase2[0] == 'K')
            {
                cout << "VLP file type detected." << endl;
                cout << "Error: Cone generator coefficients and/or duality parameters found. We only solve standard multiobjective problems. Exiting!\n";
                exit(1);
            }
        }
        fin.close();
        
        //write the problem data into either a CPLEX or glpk problem object.
        for(int i = 0; i <= max(rows,cols); i++) indices.push_back(i);
        
        cout << "The constraints:" << endl;
        for(unsigned int i = 0; i < conCoefs.size(); i++)
        {
            for(unsigned int j = 0; j < conCoefs[0].size(); j++) cout << conCoefs[i][j] << " ";
            cout << endl;
        } 
        
        cout << "The objectives:" << endl;
        for(unsigned int i = 0; i < objCoefs.size(); i++)
        {
            for(unsigned int j = 0; j < objCoefs[0].size(); j++) cout << objCoefs[i][j] << " ";
            cout << endl;
        }        
        
        #ifdef CPLEX
            char G = 'G';
            char L = 'L';
            char U = 'U';
            char B = 'B';
            char E = 'E';
            double ni = -CPX_INFBOUND;
            lp = CPXcreateprob (env,&status,argv[2]);
          	if(lp==NULL) 
          	{
            		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", 0, status);
            		exit(0);
            }
            
/*            cout << rowTypes[1] << endl;*/
            if(rowTypes[1] == 'l') status = CPXaddrows (env, lp, cols, 1, cols, &rowLB[1], &G, &intZ, indices.data(), conCoefs[0].data(), NULL, NULL);
            else if(rowTypes[1] == 'u') status = CPXaddrows (env, lp, cols, 1, cols, &rowUB[1], &L, &intZ, indices.data(), conCoefs[0].data(), NULL, NULL);
            else if(rowTypes[1] == 's') status = CPXaddrows (env, lp, cols, 1, cols, &rowLB[1], &E, &intZ, indices.data(), conCoefs[0].data(), NULL, NULL);
            for(int i = 1; i < rows; i++) 
            {
                if(rowTypes[i+1] == 'l') status = CPXaddrows (env, lp, 0, 1, cols, &rowLB[i+1], &G, &intZ, indices.data(), conCoefs[i].data(), NULL, NULL);
                else if(rowTypes[i+1] == 'u') status = CPXaddrows (env, lp, 0, 1, cols, &rowUB[i+1], &L, &intZ, indices.data(), conCoefs[i].data(), NULL, NULL);
                else if(rowTypes[i+1] == 's') status = CPXaddrows (env, lp, 0, 1, cols, &rowLB[i+1], &E, &intZ, indices.data(), conCoefs[i].data(), NULL, NULL);
            }
            
            for(int i = 0; i < cols; i++)
            {
                status = CPXchgbds (env, lp, 1, &i, &L, &ni);
                if(colTypes[i+1] == 'l') status = CPXchgbds (env, lp, 1, &i, &L, &colLB[i+1]);
                else if(colTypes[i+1] == 'u') status = CPXchgbds (env, lp, 1, &i, &U, &colUB[i+1]);
                else if(colTypes[i+1] == 's') status = CPXchgbds (env, lp, 1, &i, &B, &colLB[i+1]);
                else if(colTypes[i+1] == 'd') 
                {
                    status = CPXchgbds (env, lp, 1, &i, &L, &colLB[i+1]);
                    status = CPXchgbds (env, lp, 1, &i, &U, &colUB[i+1]);
                }
            }
            
            if(!minProb) status = CPXchgobjsen (env, lp, CPX_MAX);
            
/*            status = CPXwriteprob (env, lp, "prob_cplex.lp", "LP");*/
/*            exit(0);*/

            for(int i = 0; i < objs; i++)
            {
                lp2 = CPXcreateprob (env,&status,argv[2]);
                lp2 = CPXcloneprob (env, lp, &status);
                status = CPXchgobj(env, lp2, cols, indices.data(), objCoefs[i].data());
                myProb.AddLP(lp2);
/*                status = CPXwriteprob (env, lp, "prob_cplex.lp", "LP");*/
            }
        #else
            lp = glp_create_prob();
            
            status = glp_add_cols(lp, cols);
            status = glp_add_rows(lp, rows);
            
            for(int i = 1; i <= rows; i++)
            {
                if(rowTypes[i] == 'l') glp_set_row_bnds(lp, i, GLP_LO, rowLB[i], rowUB[i]);
                else if(rowTypes[i] == 'u') glp_set_row_bnds(lp, i, GLP_UP, rowLB[i], rowUB[i]);
                else if(rowTypes[i] == 's') glp_set_row_bnds(lp, i, GLP_FX, rowLB[i], rowUB[i]);
            }
            for(int i = 1; i <= cols; i++)
            {
                if(colTypes[i] == 'l') glp_set_col_bnds(lp, i, GLP_LO, colLB[i], colUB[i]);
                else if(colTypes[i] == 'u') glp_set_col_bnds(lp, i, GLP_UP, colLB[i], colUB[i]);
                else if(colTypes[i] == 's') glp_set_col_bnds(lp, i, GLP_FX, colLB[i], colUB[i]);
                else if(colTypes[i] == 'd') glp_set_col_bnds(lp, i, GLP_DB, colLB[i], colUB[i]);
            }
            
/*            indices.erase(indices.begin());*/
            for(int i = 0; i < rows; i++) glp_set_mat_row(lp, i+1, cols, indices.data(), conCoefs[i].data());
            
            if(minProb) glp_set_obj_dir(lp, GLP_MIN);
	        else glp_set_obj_dir(lp, GLP_MAX);
            
            for(int i = 0; i < objs; i++)
            {
                lp2 = glp_create_prob();
                glp_copy_prob(lp2, lp, GLP_ON);
                for(int j = 1; j <= cols; j++) glp_set_obj_coef(lp2, j, objCoefs[i][j]);
/*                glp_write_lp(lp2, NULL, "prob_glpk.lp");*/
/*                exit(0);*/
                myProb.AddLP(lp2);
/*                glp_write_lp(lp, NULL, "prob_glpk.lp");*/
            }
        #endif
/*        exit(0);        */
    }
    else // assume we are working with multiple *.lp files
    {
	    for(int i = 0; i < myProb.GetNumObj(); i++)
	    {
	        #ifdef CPLEX
	          	lp = CPXcreateprob (env,&status,argv[i+2]);
	          	if(lp==NULL) 
	          	{
	            		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", i+1, status);
	            		exit(0);
	            }
	            	
	            status = CPXreadcopyprob(env,lp,argv[i+2],NULL);

	          	if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(0);
	            }
	        #else
	            lp = glp_create_prob();
	        
	            status = glp_read_lp(lp, NULL, argv[i+2]);
	            if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(1);
	            }
	        #endif
	        
	        myProb.AddLP(lp);
        }
    }
    	
	/******************************************************************/
	
	/******************************************************************/
  	// Take in any other command line flags and set appropriate
  	// parameter values
  	/******************************************************************/

    if(DEBUG) cout << "setting param vals" << endl;
    myProb.SetParamVals(argc, argv);

    /******************************************************************/
	
	/******************************************************************/
	// Make sure the problems are in minimization form with less than
  	// or equal to and/or equality constraints.
  	/******************************************************************/

    myProb.SetNumRowsAndCols();
    myProb.ConvertLPs();
 
    if(DEBUG)
	{
	    #ifdef CPLEX
	        phrase = "myprob0_cplex.lp";
	    #else
		    phrase = "myprob0_glpk.lp";
		#endif
		for(int i = 0; i < myProb.GetNumObj(); i++) 
		{
			phrase[6]++;
			#ifdef CPLEX
			    status = CPXwriteprob (env, myProb.GetLP(i), phrase.c_str(), "LP");
			#else 
			    status = glp_write_lp(myProb.GetLP(i), NULL, phrase.c_str());
/*			    cout << status << endl;*/
			#endif
		}
	}

     /******************************************************************/
  

     /****************************************************
    		 Add new variables to keep track of
    		 the objective function values.
    	*****************************************************/
	
	
	if(myProb.StoreObjectivesInMainProb()) 
	{
	    if(DEBUG) cout << "adding rows for objectives" << endl;
	    myProb.AddRowsForObjectives();
	}
  	
  	/********* Initialization *****************************************/
  	
    //Implement some sort of a heuristic here that generates starting solutions.
    
/*    Simplex mySimplex(myProb.GetNumObj());*/
  				
  	/******************************************************************/

    /********* TEMPORARY CODE *****************************************/
    
    if(DEBUG) 
    {
        #ifdef CPLEX
            status = CPXwriteprob (env, myProb.GetMainLP(), "overall_prob_cplex.lp", "LP");
        #else
            status = glp_write_lp(myProb.GetMainLP(), NULL, "overall_prob_glpk.lp");
        #endif
/*        exit(0);*/
    }
  	
        if(DEBUG) cout << "finding initial simplices" << endl;
        SimplexStore t;
	myProb.DichotomicSearch(t);
  				
  	/******************************************************************/
    

  	return 0;
}


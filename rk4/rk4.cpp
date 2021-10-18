using std::vector;

// template argument dtype: data type (float, double, etc...)
// fn arguments: 
//     dfns is vector<dtype (*)(vector<dtype>)>, which is a vector of 
//     function pointers, where each fn takes a vector<dtype> as argument.
//
//     dfns F is the vector of derivative functions for each variable.
//
//     DV stores the values of the variables at each timestep, DV[0] is bc's.
//
//     h is the step size.
//
//     tmax is final time. Initial time is assumed to be 0.
template<class T>
vector<vector<T>> rk4(dfns F, vector<vector<T>> DV, T h, T tmax){
    int numSteps = (int)((tmax - DV[0][0])/h);

    //do rk4
    for(int i = 1; i <= numSteps+1; i++){
        vector<T> next; 
        vector<T> arg;
        vector<vector<T>> K = {{0},{0},{0},{0}};
        //fill next values with 0's
        for(unsigned j = 0; j < DV[0].size()-1; j++){next.push_back(0);}
        for(unsigned j = 0; j < 4; j++){K[j] = next;}
        next.push_back(0);
        next[0] = i*h; //next time value = i*stepsize
        for(unsigned k = 0; k < 4; k++){
            arg = DV[i-1];
            // make rk4 arguments to derivative functions
            if(k == 3){
                arg[0] += h;
                for(unsigned n = 1; n <= K[0].size(); n++){
                    arg[n] += K[2][n-1];
                }
            }
            else if(k > 0){
                arg[0] += (h/2);
                for(unsigned n = 1; n <= K[0].size(); n++){
                    arg[n] += K[k-1][n-1]/2;
                }
            }
            //call derivative fn to get rk4 coefficient
            for(unsigned j = 0; j < F.size(); j++){
                K[k][j] = h*F[j](arg);
            }
        }
        for(unsigned j = 1; j < DV[0].size(); j++){
            next[j]=DV[i-1][j]+(K[0][j-1]+2*K[1][j-1]+2*K[2][j-1]+K[3][j-1])/6;
        }
        DV.push_back(next);
    }
    return DV;
}

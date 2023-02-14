#include "interpolate.h"

void interpolate(const real &tNext, const std::vector<real> &tVecForInterp, const std::vector< std::vector<real> > &xIntegForInterp, Simulation &sim){
    size_t tLen = tVecForInterp.size();
    size_t numStates = 6*sim.integParams.nInteg;
    std::vector< std::vector< std::vector<real> > > c(numStates, std::vector< std::vector<real> >(tLen, std::vector<real>(tLen, 0.0)));
    std::vector< std::vector<real> > coeffs(numStates, std::vector<real>(tLen, 0.0));

    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            c[i][j][0] = xIntegForInterp[j][i];
        }
        for (size_t j = 1; j < tLen; j++){
            for (size_t k = 0; k < tLen-j; k++){
                c[i][k][j] = (c[i][k+1][j-1] - c[i][k][j-1])/(tVecForInterp[k+j]-tVecForInterp[k]);
            }
        }
    }

    for (size_t i = 0; i < numStates; i++){
        for (size_t j = 0; j < tLen; j++){
            coeffs[i][j] = c[i][0][j];
        }
    }

    static size_t interpIdx = 0;
    if (interpIdx == sim.tEval.size()){
        interpIdx = 0;
    }
    bool forwardIntegrate = tVecForInterp[0] < tVecForInterp[tLen-1];
    bool backwardIntegrate = tVecForInterp[0] > tVecForInterp[tLen-1];
    while ( interpIdx < sim.tEval.size()
            &&( (forwardIntegrate && (sim.tEval[interpIdx] == tVecForInterp[0] || (sim.tEval[interpIdx] > tVecForInterp[0] && sim.tEval[interpIdx] <= tNext)))
            ||  (forwardIntegrate && sim.tEval[interpIdx] <= sim.integParams.t0 && sim.tEval[interpIdx]+sim.tEvalMargin >= sim.integParams.t0) || (forwardIntegrate && sim.tEval[interpIdx] >= sim.integParams.tf && sim.tEval[interpIdx]-sim.tEvalMargin <= sim.integParams.tf)
            ||  (backwardIntegrate && (sim.tEval[interpIdx] == tVecForInterp[0] || (sim.tEval[interpIdx] < tVecForInterp[0] && sim.tEval[interpIdx] >= tNext)))
            ||  (backwardIntegrate && sim.tEval[interpIdx] >= sim.integParams.t0 && sim.tEval[interpIdx]-sim.tEvalMargin <= sim.integParams.t0) || (backwardIntegrate && sim.tEval[interpIdx] <= sim.integParams.tf && sim.tEval[interpIdx]+sim.tEvalMargin >= sim.integParams.tf) )
        ){
        real tInterp;
        for (size_t i = 0; i < tLen; i++){
            tInterp = sim.tEval[interpIdx];
            if (tInterp == tVecForInterp[i]){
                sim.xIntegEval.push_back(xIntegForInterp[i]);
                interpIdx++;
                // std::cout << "exactly interpolated tInterp = " << tInterp << std::endl;
                continue;
            }
        }
        tInterp = sim.tEval[interpIdx];
        // std::cout << "tInterp = " << tInterp << std::endl;
        std::vector<real> xInterp(numStates, 0.0);
        size_t n = tLen-1;
        for (size_t i = 0; i < numStates; i++){
            xInterp[i] = coeffs[i][n];
        }
        for (size_t i = 1; i < n+1; i++){
            for (size_t j = 0; j < numStates; j++){
                xInterp[j] = coeffs[j][n-i] + (tInterp-tVecForInterp[n-i])*xInterp[j];
            }
        }
        sim.xIntegEval.push_back(xInterp);
        interpIdx++;
    }
}

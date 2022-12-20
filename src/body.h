#ifndef BODY_H
#define BODY_H

typedef class Body{
    private:

    public:
    real t0;
    real mass;
    real radius;
    std::string name;
    std::vector<real> pos;
    std::vector<real> vel;
    bool isNongrav=false;
} body;

typedef class SpiceBody: public Body{
    private:

    public:
    int spiceID;
    const bool isSpice=true;
    // constructor
    SpiceBody(std::string name, int spiceID, real t0, real mass, real radius);

} spiceBody;

typedef class IntegBody: public Body{
    private:

    public:
    const bool isInteg=true;
    std::vector< std::vector<real> > covariance;
    NongravParams ngParams;
    // constructor
    IntegBody(std::string name, real t0, real mass, real radius, std::vector<real> pos, std::vector<real> vel, std::vector< std::vector<real> > covariance, real a1, real a2, real a3, real alpha, real k, real m, real n, real r0_au);

} integBody;

#endif

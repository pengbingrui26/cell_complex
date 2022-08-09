#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Eigen>

const double tol = 1e-4;


// class ==============================================================

class lattvec
{
public:
    Eigen::MatrixXd array;
    double get_norm(void);
    bool is_zero(void);
    bool operator==(lattvec v1);
    lattvec home(Eigen::MatrixXd AA, Eigen::MatrixXd BB);
    int order(Eigen::MatrixXd AA, Eigen::MatrixXd BB);
};

double lattvec::get_norm(void)
{
    return this->array.norm();
}

bool lattvec::is_zero(void)
{
    double maxcoeff = std::abs(this->array.maxCoeff());
    double mincoeff = std::abs(this->array.minCoeff());
    double max_abs_element = std::max(maxcoeff, mincoeff);
    return max_abs_element < tol;
}

bool lattvec::operator==(lattvec v1)
{
    lattvec tmp;
    tmp.array = this->array - v1.array;
    return tmp.is_zero();
}

lattvec lattvec::home(Eigen::MatrixXd AA, Eigen::MatrixXd BB)
{
    // Move the vector to home cell
    Eigen::MatrixXd direc = this->array * BB.transpose();

    if (this->array.rows() == 1)
    {
        for (int ii = 0; ii < direc.cols(); ii++)
        {
            if (std::abs(direc(0, ii) - std::round(direc(0, ii))) < tol)
            {
                direc(0, ii) = 0.;
            }
            else
            {
                direc(0, ii) -= std::floor(direc(0, ii));
            }
        }
    }
    else if (this->array.rows() == 2)
    {
        for (int ii = 0; ii < direc.rows(); ii++)
        {
            for (int jj = 0; jj < direc.cols(); jj++)
            {
                if (std::abs(direc(ii, jj) - std::round(direc(ii, jj))) < tol)
                {
                    direc(ii, jj) = 0.;
                }
                else
                {
                    direc(ii, jj) -= std::floor(direc(ii, jj));
                }
            }
        }
    }
    lattvec vv;
    vv.array = direc;
    return vv;
}

int lattvec::order(Eigen::MatrixXd AA, Eigen::MatrixXd BB)
{
    int maxorder = 10;
    bool find_order = 0;
    Eigen::MatrixXd direc = this->array * BB.transpose();

    int ii = 1;
    for (; ii < maxorder + 1; ii++)
    {
        lattvec tmp;
        tmp.array = ii * direc;
        if (tmp.home(AA, BB).is_zero())
        {
            find_order = 1;
            break;
        }
    }

    if (!find_order)
    {
        throw "order not found !!!";
    }
    return ii;
}

// ===============================================================================================

class symop
{
public:
    lattvec mat;
    lattvec tau;
    lattvec matP;
    lattvec tauP;

    double det;
    double alpha;
    int nfold;
    lattvec axis;

    Eigen::MatrixXd AA;
    Eigen::MatrixXd BB;
    Eigen::MatrixXd AAC;
    Eigen::MatrixXd BBC;

    symop(lattvec, lattvec, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);
};

symop::symop(lattvec mat, lattvec tau, Eigen::MatrixXd AA, Eigen::MatrixXd BB, Eigen::MatrixXd AAC, Eigen::MatrixXd BBC)
{
    this->mat = mat;
    this->tau = tau;
    this->AA = AA;
    this->BB = BB;
    this->AAC = AAC;
    this->BBC = BBC;

    // check cell
    Eigen::MatrixXd iden_matr;
    iden_matr << 1., 0., 0., 
                 0., 1., 0.,
                 0., 0., 1.;

    Eigen::MatrixXd BBAA = this->BB.transpose() * this->AA - iden_matr;
    double maxabscoeff_BBAA = BBAA.array().abs().maxCoeff();

    Eigen::MatrixXd BBAAC = this->BBC.transpose() * this->AAC - iden_matr;
    double maxabscoeff_BBAAC = BBAA.array().abs().maxCoeff();

    if (maxabscoeff_BBAA > 1e-4 || maxabscoeff_BBAAC > 1e-4)  
    {
        throw "wrong reciprocal lattice !!!";        
    }

    //transform matrix and translation vector
    this->tau = tau.home(AA, BB);
    this->matP.array = this->BB * this->mat.array * this->AA.transpose();
    this->tauP.array = this->BB * this->tau.array;
    //this->matP.array = std::round(this->matP);
}
;



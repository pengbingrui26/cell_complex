#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include "spg.cpp"


// test =========================================================

void test_array()
{ 
    Eigen::MatrixXd H_tmp(4, 4); 
    double t = 1.;  
    double U = 2.;  
    H_tmp << 0., 0., -t, -t,
            0., 0., -t, -t,
            -t, -t, U, 0.,
            -t, -t, 0., U;  
    lattvec hh;  
    hh.array = H_tmp;  
    std::cout << hh.array << std::endl;  
    std::cout << hh.get_norm() << std::endl;

    Eigen::MatrixXd VV(1, 4);
    VV << 0., 0., -t, -t;
    lattvec vv;
    vv.array = VV;
    std::cout << vv.array << std::endl;
    std::cout << vv.get_norm() << std::endl;
}

void test_iszero()
{
    Eigen::MatrixXd H_tmp(4, 4);
    double t = 0.;
    double U = 0.;
    H_tmp << 0., 0., -t, -t,
        0., 0., -t, -t,
        -t, -t, U, 0.,
        -t, -t, 0., U;
    lattvec hh;
    hh.array = H_tmp;
    std::cout << hh.is_zero() << std::endl;
}

void test_equal()
{
    Eigen::MatrixXd v1(1, 3);
    v1 << 1., 2., 3.;
    Eigen::MatrixXd v2(1, 3);
    v2 << 1., 4., 3.;
    lattvec arr1;
    arr1.array = v1;
    lattvec arr2;
    arr2.array = v2;
    std::cout << (arr1 == arr2) << std::endl;
}

void test_home()
{
    Eigen::MatrixXd standard_AA(3, 3);
    Eigen::MatrixXd standard_BB(3, 3);

    standard_AA << 1., 0., 0.,
        0., 1., 0.,
        0., 0., 1;

    standard_BB << 1., 0., 0.,
        0., 1., 0.,
        0., 0., 1;

    Eigen::MatrixXd vv(1, 3);
    vv << 1, 2.9, 3;
    lattvec arr;
    arr.array = vv;

    lattvec arr_new = arr.home(standard_AA, standard_BB);
    std::cout << arr_new.array << std::endl;
}

void test_order()
{
    Eigen::MatrixXd standard_AA(3, 3);
    Eigen::MatrixXd standard_BB(3, 3);

    standard_AA << 1., 0., 0.,
        0., 1., 0.,
        0., 0., 1;

    standard_BB << 1., 0., 0.,
        0., 1., 0.,
        0., 0., 1;

    Eigen::MatrixXd vv(1, 3);
    vv << 0, 0, 0.25;
    lattvec arr;
    arr.array = vv;

    int n = arr.order(standard_AA, standard_BB);
    std::cout << n << std::endl;
}

// ==============================================================

int main()
{
    // test_array();
    // test_iszero();
    // test_equal();
    // test_home();
    test_order();
}


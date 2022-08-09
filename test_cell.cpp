# include <iostream>
# include <vector>
# include <numeric>
# include <cmath>
# include <algorithm>
# include <assert.h>
# include <tuple>
# include "cell.cpp"


// test ===================================================================

void test_print()
{
    std::vector<int> vv = {1,1,1,3};
    MyPrint(vv);

}


void test_slice()
{
    std::vector<int> vv = {1,2,3,4,5,6};
    MyPrint(Slice(vv, 0, 2));
}


void test_get_inner_product()
{
    std::vector<int> vv = {1,2,3};
    std::cout << get_inner_product(vv, vv) << std::endl;
}


void test_vector_add()
{
    std::vector<int> v1 = {1,2,3};
    std::vector<int> v2 = {1,2,3};
    std::vector<int> v3 = vector_add(v1, v2);
    MyPrint(v3);
}


void test_vector_equal()
{
    std::vector<int> v1 = {1,2,3};
    std::vector<int> v2 = {1,2,3,5};
    std::cout << vector_equal(v1, v2) << std::endl; 
}


void test_vector_in()
{
    std::vector<int> v = {1,2,3};
    std::vector<std::vector<int>> vv = {{1,2,3}, {4,5,6}};
    std::cout << vector_in(vv, v) << std::endl; 
}

void test_get_segment()
{
    std::vector<std::vector<double>> ineq2 = {{1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}};
    auto seg = get_segment(ineq2, 0);
    //std::vector<double> xmin, xmax;
    auto v1 = std::get<0>(seg);
    auto v2= std::get<1>(seg);
    //std::cout << xmin << ", " << xmax << std::endl;
    MyPrint(v1);
    MyPrint(v2);
}


// ====================================================================

int main()
{
    //test_print();
    //test_slice();
    //test_get_inner_product();
    //test_vector_add();
    //test_vector_equal();
    test_vector_in();
    //test_get_segment();
}



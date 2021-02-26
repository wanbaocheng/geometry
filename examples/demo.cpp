#include <iostream>
#include "geometry.hpp"

using namespace std;

int main() {
    geometry::Rigid3<double> r3_1;
    Eigen::AngleAxisd aa(r3_1.rotation());
    cout << r3_1.rotation()(0, 0) << endl;
    cout << Eigen::AngleAxisd(r3_1.rotation()).angle() << endl;

    geometry::Rigid3<float> r3_2(geometry::Vector3<float>{3.14159, -1.42233, 1.59033},
                                 geometry::EulerConvention::xyz_s);
    cout << r3_2.inverse().eulerAngles(geometry::EulerConvention::xyz_s) << endl;
}
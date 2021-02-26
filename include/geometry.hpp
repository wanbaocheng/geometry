//
// Author: Wan Baocheng
// Email: wanbaocheng1@jd.com
// Date:  2020/01/16
//
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <Eigen/Geometry>

namespace geometry {
    constexpr double PI = M_PI;

    template<template<typename> class DerivedType, typename Scalar>
    class Angle {
    public:
        Angle() = default;

        explicit Angle(const Scalar &value) {
            this->set(value);
        }

        Scalar get() const {
            return this->value_;
        }

        void set(const Scalar &val) {
            this->value_ = val;
        }


        Angle(const Angle &angle) {
            this->set(angle.get());
        }

        DerivedType<Scalar> operator-() const {
            return DerivedType<Scalar>{-value_};
        }

        DerivedType<Scalar> operator+(const DerivedType<Scalar> &other) const {
            return DerivedType<Scalar>{value_ + other.value_};
        }

        DerivedType<Scalar> operator-(const DerivedType<Scalar> &other) const {
            return DerivedType<Scalar>{value_ - other.value_};
        }

        DerivedType<Scalar> &operator+=(const DerivedType<Scalar> &other) {
            this->set(value_ + other.value_);
            return static_cast<DerivedType<Scalar> &>(*this);
        }

        DerivedType<Scalar> &operator-=(const DerivedType<Scalar> &other) {
            this->set(value_ - other.value_);
            return static_cast<DerivedType<Scalar> &>(*this);
        }

        DerivedType<Scalar> operator*(const Scalar scalar) const {
            return DerivedType<Scalar>{value_ * scalar};
        }

        DerivedType<Scalar> operator/(const Scalar scalar) const {
            return DerivedType<Scalar>{value_ / scalar};
        }

        DerivedType<Scalar> &operator*=(const Scalar scalar) {
            this->set(value_ * scalar);
            return static_cast<DerivedType<Scalar> &>(*this);
        }

        DerivedType<Scalar> &operator/=(const Scalar scalar) {
            this->set(value_ / scalar);
            return static_cast<DerivedType<Scalar> &>(*this);
        }

        bool operator==(const DerivedType<Scalar> &other) const {
            return value_ == other.value_;
        }

        bool operator!=(const DerivedType<Scalar> &other) const {
            return value_ != other.value_;
        }

        bool operator<=(const DerivedType<Scalar> &other) const {
            return value_ <= other.value_;
        }

        bool operator>=(const DerivedType<Scalar> &other) const {
            return value_ >= other.value_;
        }

        bool operator<(const DerivedType<Scalar> &other) const {
            return value_ < other.value_;
        }

        bool operator>(const DerivedType<Scalar> &other) const {
            return value_ > other.value_;
        }

    private:
        Scalar value_ = Scalar{0};
    };

    template<typename Scalar>
    class Deg;

    template<typename Scalar>
    class Rad final : public Angle<Rad, Scalar> {
    public:
        using Base = Angle<Rad, Scalar>;
        using Base::Base;

        Rad(const Deg<Scalar> &deg) {
            this->set(deg.get() * PI / 180);
        }

        friend std::ostream &operator<<(std::ostream &out, const Rad<Scalar> &rad) {
            out << rad.get();
            return out;
        }
    };

    template<typename Scalar>
    class Deg final : public Angle<Deg, Scalar> {
    public:
        using Base = Angle<Deg, Scalar>;
        using Base::Base;

        Deg(const Rad<Scalar> &rad) {
            this->set(rad.get() * 180 / PI);
        }

        friend std::ostream &operator<<(std::ostream &out, const Deg<Scalar> &deg) {
            out << deg.get() << "°";
            return out;
        }
    };

    template<typename _Scalar, int _Dim>
    using Vector = Eigen::Matrix<_Scalar, _Dim, 1>;
    template<typename _Scalar>
    using Vector2 = Vector<_Scalar, 2>;
    using Vector2f = Vector2<float>;
    using Vector2d = Vector2<double>;
    template<typename _Scalar>
    using Vector3 = Vector<_Scalar, 3>;
    using Vector3f = Vector3<float>;
    using Vector3d = Vector3<double>;
    template<typename _Scalar>
    using Vector4 = Vector<_Scalar, 4>;
    using Vector4f = Vector4<float>;
    using Vector4d = Vector4<double>;

    template<typename _Scalar>
    class Rigid2 {
        using Scalar = _Scalar;
    public:
        Rigid2() {
            this->transform = Eigen::Transform<_Scalar, 2, Eigen::Affine>::Identity();
        }

        Rigid2(const Rigid2 &rigid) {
            this->transform = rigid.transform;
        }

        explicit Rigid2(const Scalar &angle_rad) {
            this->transform = Eigen::Rotation2D<Scalar>(angle_rad);
        }

        explicit Rigid2(const Vector2<Scalar> &translate) {
            this->transform = Eigen::Translation<Scalar, 2>(translate[0], translate[1]);
        }

        Rigid2(const Scalar &angle_rad, const Vector2<Scalar> &translate) {
            this->transform = Eigen::Translation<Scalar, 2>(translate[0], translate[1]);
            this->transform *= Eigen::Rotation2D<Scalar>(angle_rad);
        }

        Eigen::Matrix<Scalar, 2, 2> rotation() const {
            return this->transform.rotation();
        }

        Rigid2<Scalar> &rotate(const Scalar &angle_rad) {
            this->transform *= Eigen::Rotation2D<Scalar>(angle_rad);
            return *this;
        }

        Scalar angle() const {
            Eigen::Rotation2D<Scalar> r2d;
            r2d = this->transform.rotation();
            return r2d.angle();
        }

        Vector2<Scalar> translation() const {
            return this->transform.translation();
        }

        Rigid2<Scalar> &translate(const Vector2<Scalar> &translate) {
            this->transform.translate(translate);
            return *this;
        }

        Rigid2<Scalar> &operator*=(const Rigid2<Scalar> &rigid) {
            this->transform = this->transform * rigid.transform;
            return *this;
        }

        Rigid2<Scalar> operator*(const Rigid2<Scalar> &rigid) {
            Rigid2<Scalar> rigid2{};
            rigid2.transform = this->transform * rigid.transform;
            return rigid2;
        }

        Rigid2 inverse() const {
            Rigid2 rigid{*this};
            rigid.transform = rigid.transform.inverse();
            return rigid;
        }

        static Rigid2 Identity() {
            return Rigid2{};
        }

    private:
        Eigen::Transform<_Scalar, 2, Eigen::Affine> transform;
    };

    enum class EulerConvention : unsigned int {
        xyx_r = 0b00010000, xzx_r = 0b00100000, yxy_r = 0b01000100, yzy_r = 0b01100100,
        zxz_r = 0b10001000, zyz_r = 0b10011000, xyz_r = 0b00011000, xzy_r = 0b00100100,
        yxz_r = 0b01001000, yzx_r = 0b01100000, zxy_r = 0b10000100, zyx_r = 0b10010000,  // rotation axis convention
        xyx_s = 0b00010001, xzx_s = 0b00100001, yxy_s = 0b01000101, yzy_s = 0b01100101,
        zxz_s = 0b10001001, zyz_s = 0b10011001, xyz_s = 0b00011001, xzy_s = 0b00100101,
        yxz_s = 0b01001001, yzx_s = 0b01100001, zxy_s = 0b10000101, zyx_s = 0b10010001   // static axis convention
    };

    template<typename _Scalar>
    class Rigid3 {
    public:
        using Scalar = _Scalar;

        Rigid3() {
            this->transform = Eigen::Transform<_Scalar, 3, Eigen::Affine>::Identity();
        }

        Rigid3(const Rigid3 &rigid) {
            this->transform = rigid.transform;
        }

        Rigid3(const Vector3<Scalar> &eulerAngles, EulerConvention order) {
            this->transform = Rigid3::eulerAngles2Transform(eulerAngles, order);
        }

        Rigid3(const Scalar &angle, Vector3<Scalar> axis) {
            this->transform = Eigen::AngleAxis<Scalar>(angle, axis);
        }

        Rigid3(const Deg<Scalar> &deg, Vector3<Scalar> axis) {
            this->transform = Eigen::AngleAxis<Scalar>(RadD(deg).get(), axis);
        }

        Rigid3(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) {
            this->transform = Eigen::Quaternion<Scalar>(w, x, y, z);
        }

        Rigid3(const Vector4<Scalar> &xyzw) {
            this->transform = Eigen::Quaternion<Scalar>(xyzw);
        }

        explicit Rigid3(const Vector3<Scalar> &translate) {
            this->transform = Eigen::Translation<Scalar, 3>(translate[0], translate[1], translate[2]);
        }

        Rigid3(const Vector3<Scalar> &eulerAngles, EulerConvention order, const Vector3<Scalar> &translate) {
            this->transform = Eigen::Translation<Scalar, 3>(translate[0], translate[1], translate[2]);
            this->transform *= Rigid3::eulerAngles(eulerAngles, order);
        }

        Vector3<Scalar> eulerAngles(EulerConvention order) const {
            auto order_num = (unsigned int) order;
            if ((order_num & 0b11u) == 0) {
                return this->transform.rotation().eulerAngles((order_num >> 6u) & 0b11u,
                                                              (order_num >> 4u) & 0b11u,
                                                              (order_num >> 3u) & 0b11u);
            } else {
                return this->transform.rotation().eulerAngles((order_num >> 2u) & 0b11u,
                                                              (order_num >> 4u) & 0b11u,
                                                              (order_num >> 6u) & 0b11u);
            }
        }

        Vector4<Scalar> quaternion() const {
            Eigen::Quaternion<Scalar> q = this->transform;
            return q.coeffs();
        }

        Eigen::Matrix<Scalar, 3, 3> rotation() const {
            return this->transform.rotation();
        }

        Rigid3<Scalar> &rotate(const Vector4<Scalar> &xyzw) {
            this->transform *= Eigen::Quaternion<Scalar>(xyzw);
            return *this;
        }

        Rigid3<Scalar> &rotate(const Scalar &w, const Scalar &x, const Scalar &y, const Scalar &z) {
            this->transform = *Eigen::Quaternion<Scalar>(w, x, y, z);
            return *this;
        }

        Vector3<Scalar> translation() const {
            return this->transform.translation();
        }

        void translation(const Vector3<Scalar> &trans) {
            this->transform.translation() = trans;
        }

        Rigid3<Scalar> &translate(const Vector3<Scalar> &translate) {
            this->transform.translate(translate);
            return *this;
        }

        Rigid3<Scalar> &pretranslate(const Vector3<Scalar> &translate) {
            this->transform.pretranslate(translate);
            return *this;
        }

        Rigid3<Scalar> &setTranslate(const Vector3<Scalar> &translate) {
            this->transform.translate(translate);
            return *this;
        }

        Rigid3<Scalar> &operator*=(const Rigid3<Scalar> &rigid) {
            this->transform = this->transform * rigid.transform;
            return *this;
        }

        Rigid3<Scalar> operator*(const Rigid3<Scalar> &rigid) {
            Rigid3<Scalar> rigid3{};
            rigid3.transform = this->transform * rigid.transform;
            return rigid3;
        }

        Rigid3 inverse() const {
            Rigid3 rigid{*this};
            rigid.transform = rigid.transform.inverse();
            return rigid;
        }

        static Rigid3 Identity() {
            return Rigid3{};
        }

        static const Vector3<Scalar> X;
        static const Vector3<Scalar> Y;
        static const Vector3<Scalar> Z;

    private:
        Eigen::Transform<_Scalar, 3, Eigen::Affine> transform;

        static Eigen::Transform<_Scalar, 3, Eigen::Affine>
        eulerAngles2Transform(const Vector3<Scalar> &eulerAngles, EulerConvention order) {
            Eigen::Transform<_Scalar, 3, Eigen::Affine> transform;
            auto order_num = (unsigned int) order;
            if ((order_num & 0b11u) == 0) {
                transform = Eigen::AngleAxis<Scalar>(eulerAngles[0], Vector3<Scalar>::Unit((order_num >> 6u) & 0b11u)) *
                            Eigen::AngleAxis<Scalar>(eulerAngles[1], Vector3<Scalar>::Unit((order_num >> 4u) & 0b11u)) *
                            Eigen::AngleAxis<Scalar>(eulerAngles[2], Vector3<Scalar>::Unit((order_num >> 2u) & 0b11u));
            } else {
                transform = Eigen::AngleAxis<Scalar>(eulerAngles[0], Vector3<Scalar>::Unit((order_num >> 2u) & 0b11u)) *
                            Eigen::AngleAxis<Scalar>(eulerAngles[1], Vector3<Scalar>::Unit((order_num >> 4u) & 0b11u)) *
                            Eigen::AngleAxis<Scalar>(eulerAngles[2], Vector3<Scalar>::Unit((order_num >> 6u) & 0b11u));
            }
            return transform;
        }
    };

    template<typename _Scalar>
    const Vector3<_Scalar> Rigid3<_Scalar>::X = Eigen::Matrix<_Scalar,
            3, 1>::UnitX();
    template<typename _Scalar>
    const Vector3<_Scalar> Rigid3<_Scalar>::Y = Eigen::Matrix<_Scalar,
            3, 1>::UnitY();
    template<typename _Scalar>
    const Vector3<_Scalar> Rigid3<_Scalar>::Z = Eigen::Matrix<_Scalar,
            3, 1>::UnitZ();

    template<typename _Scalar>
    class Transform3 {
    public:
        using Scalar = _Scalar;

        Transform3(std::string id_frame_source,
                   std::string id_frame_target,
                   const Rigid3<Scalar> &rigid3) :
                id_frame_source(std::move(id_frame_source)),
                id_frame_target(std::move(id_frame_target)),
                rigid3(rigid3) {}

        Transform3(const Transform3(&transform)) :
                id_frame_source(transform.id_frame_source),
                id_frame_target(transform.id_frame_target),
                rigid3(transform.rigid3) {}

        const std::string id_frame_source;
        const std::string id_frame_target;
        Rigid3<_Scalar> rigid3;
    };

    template<typename _Scalar, typename TimeStamped>
    class Transform3Stamped {
        using Scalar = _Scalar;
    public:
        Transform3Stamped() = default;

        Transform3Stamped(std::string id_frame_source,
                          std::string id_frame_target,
                          const Rigid3<Scalar> &rigid3,
                          TimeStamped time_stamped) :
                id_frame_source(std::move(id_frame_source)),
                id_frame_target(std::move(id_frame_target)),
                rigid3(rigid3),
                time_stamped(time_stamped) {}

        Transform3Stamped(const Transform3Stamped &transform_stamped) :
                id_frame_source(transform_stamped.id_frame_source),
                id_frame_target(
                        transform_stamped.id_frame_target),
                rigid3(transform_stamped.rigid3),
                time_stamped(transform_stamped.time_stamped) {}

        Transform3Stamped &operator=(const Transform3Stamped &transform_stamped) {
            if (this == &transform_stamped)
                return *this;
            this->id_frame_target = transform_stamped.id_frame_target;
            this->id_frame_source = transform_stamped.id_frame_source;
            this->rigid3 = transform_stamped.rigid3;
            this->time_stamped = transform_stamped.time_stamped;
        }

        std::string id_frame_source;
        std::string id_frame_target;
        TimeStamped time_stamped;
        Rigid3<_Scalar> rigid3;
    };
}

#endif

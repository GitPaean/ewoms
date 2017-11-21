/*
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_POINT2D_HEADER_INCLUDED
#define OPM_POINT2D_HEADER_INCLUDED

#include <array>

namespace Opm {



    namespace detail {


        template <class ValueType>
        class Point2D
        {
            public:
                Point2D(const ValueType& xi, const ValueType& yi)
                    : x_(xi),
                      y_(yi) {
                      }

                Point2D()
                    : x_(0.),
                      y_(0.) {
                      }

                const ValueType& getX() const
                {
                    return x_;
                }

                const ValueType& getY() const
                {
                    return y_;
                }

                void setX(const ValueType& x)
                {
                    x_ = x;
                }

                void setY(const ValueType& y)
                {
                    y_ = y;
                }

                /// Finding the intersection point of a line segment and a line.
                /// return true, if found.
                static bool findIntersection(const std::array<Point2D<ValueType>, 2> line_segment1,
                                             const std::array<Point2D<ValueType>, 2> line2, Point2D<ValueType>& intersection_point)
                {

                    const ValueType& x1 = line_segment1[0].getX();
                    const ValueType& y1 = line_segment1[0].getY();
                    const ValueType& x2 = line_segment1[1].getX();
                    const ValueType& y2 = line_segment1[1].getY();

                    const ValueType& x3 = line2[0].getX();
                    const ValueType& y3 = line2[0].getY();
                    const ValueType& x4 = line2[1].getX();
                    const ValueType& y4 = line2[1].getY();

                    const ValueType& d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

                    if (d == 0.) {
                        return false;
                    }

                    const ValueType& x = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
                    const ValueType& y = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

                    if (x >= Opm::min(x1, x2) && x <= Opm::max(x1, x2)) {
                        intersection_point.setX(x);
                        intersection_point.setY(y);
                        return true;
                    } else {
                       return false;
                    }
                }

            private:
                ValueType x_;
                ValueType y_;

        }; // class Point2D


    } // namespace detail



} // namespace Opm

#endif // OPM_POINT2D_HEADER_INCLUDED


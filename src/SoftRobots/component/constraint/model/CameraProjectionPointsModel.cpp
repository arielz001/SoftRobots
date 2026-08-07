/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
*                           Plugin SoftRobots                                 *
*                                                                             *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
******************************************************************************/
#define SOFTROBOTS_CAMERAPROJECTIONMODEL_CPP
#include <SoftRobots/component/constraint/model/CameraProjectionPointsModel.inl>
#include <sofa/core/ObjectFactory.h>

namespace softrobots::constraint
{

using namespace sofa::defaulttype;
using sofa::core::ConstraintParams;

template<> SOFA_SOFTROBOTS_API
void CameraProjectionPointsModel<Rigid3Types>::normalizeDirections()
{
    VecDeriv directions;
    directions.resize(6);
    for(unsigned int i = 0; i < 6; i++)
    {
        directions[i] = d_directions.getValue()[i];
        Vec<3, Real> vector1 {directions[i][0], directions[i][1], directions[i][2]};
        Vec<3, Real> vector2 {directions[i][3], directions[i][4], directions[i][5]};
        vector1.normalize();
        vector2.normalize();
        directions[i] = Deriv(vector1, vector2);
    }
    d_directions.setValue(directions);
}

int CameraProjectionPointsModelClass = sofa::core::RegisterObject("This component computes the 2D pinhole projection (u, v) of 3D point(s) onto a camera image plane.")
                .add< CameraProjectionPointsModel<Vec1Types> >()
                .add< CameraProjectionPointsModel<Vec2Types> >()
                .add< CameraProjectionPointsModel<Vec3Types> >()
                .add< CameraProjectionPointsModel<Rigid3Types> >(true);

template<> SOFA_SOFTROBOTS_API
void CameraProjectionPointsModel<Vec1Types>::drawPoints(const VisualParams* vparams, const std::vector<Coord> &points, float size, const RGBAColor& color) {
    SOFA_UNUSED(vparams);
    SOFA_UNUSED(points);
    SOFA_UNUSED(size);
    SOFA_UNUSED(color);
}

template<> SOFA_SOFTROBOTS_API
void CameraProjectionPointsModel<Vec3Types>::drawPoints(const VisualParams* vparams, const std::vector<Coord> &points, float size, const RGBAColor& color) {
    vparams->drawTool()->drawPoints(points, size, color);
}

template<> SOFA_SOFTROBOTS_API
void CameraProjectionPointsModel<Vec2Types>::drawPoints(const VisualParams* vparams, const std::vector<Coord> &points, float size, const RGBAColor& color) {
    sofa::type::vector<Vec3> pointsVec3;
    pointsVec3.reserve(points.size());
    for (const auto& point : points)
        pointsVec3.push_back(Vec3(point[0], point[1], 0.0));
    vparams->drawTool()->drawPoints(pointsVec3, size, color);
}

template<> SOFA_SOFTROBOTS_API
void CameraProjectionPointsModel<Rigid3Types>::drawPoints(const VisualParams* vparams, const std::vector<Coord> &points, float size, const RGBAColor& color) {
    sofa::type::vector<Vec3> pointsVec3;
    pointsVec3.reserve(points.size());
    for (const auto& point : points)
        pointsVec3.push_back(point.getCenter());
    vparams->drawTool()->drawPoints(pointsVec3, size, color);
}

template class SOFA_SOFTROBOTS_API CameraProjectionPointsModel<Vec1Types>;
template class SOFA_SOFTROBOTS_API CameraProjectionPointsModel<Vec2Types>;
template class SOFA_SOFTROBOTS_API CameraProjectionPointsModel<Vec3Types>;
template class SOFA_SOFTROBOTS_API CameraProjectionPointsModel<Rigid3Types>;

} // namespace softrobots::constraint
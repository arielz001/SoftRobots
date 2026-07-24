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
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
******************************************************************************/
#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <SoftRobots/component/constraint/model/CameraProjectionModel.h>
#include <Eigen/Dense>
#include <array>

using namespace Eigen;

namespace softrobots::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;
using sofa::core::ConstVecCoordId;
using sofa::helper::WriteAccessor;
using sofa::helper::ReadAccessor;
using sofa::type::vector;
using sofa::type::Vec;
using sofa::type::Vec3;
using sofa::type::RGBAColor;

// Constructor
template<class DataTypes>
CameraProjectionModel<DataTypes>::CameraProjectionModel(MechanicalState* object)
    : Inherit1(object)
    , d_indices(initData(&d_indices, "indices", "Indices of nodes to project"))
    , d_focalLength(initData(&d_focalLength, sofa::type::Vec2(600.0, 600.0), "focalLength", "Focal length [fx, fy] in pixels."))
    , d_principalPoint(initData(&d_principalPoint, sofa::type::Vec2(320.0, 240.0), "principalPoint", "Principal point [u0, v0] in pixels."))
    , d_ellipseRadius(initData(&d_ellipseRadius, 0.05, "ellipseRadius", "Real 3D radius of the target ellipse in meters."))
    , d_weight(initData(&d_weight, sofa::type::vector<Real>(Deriv::total_size, 1.), "weight", "Weight for minimization"))
    , d_directions(initData(&d_directions, "directions", "Directions vector"))
    , d_Jacobian(initData(&d_Jacobian, "Jacobian", "Jacobian relating node motion to 2D ellipse changes"))
    , d_useDirections(initData(&d_useDirections, "useDirections", "Directions selection mask"))
    , d_delta(initData(&d_delta, "delta", "Distance to target"))
{
    d_delta.setReadOnly(true);

    this->addUpdateCallback("updateWeight", {&d_weight}, [this](const sofa::core::DataTracker& t)
    {
        SOFA_UNUSED(t);
        auto weight = sofa::helper::getWriteAccessor(d_weight);
        if (weight.size() != Deriv::total_size)
        {
            Real w = weight.empty() ? 1. : weight[0];
            d_weight.setValue(sofa::type::vector<Real>(Deriv::total_size, w));
        }
        return sofa::core::objectmodel::ComponentState::Valid;
    }, {});
}

// Destructor
template<class DataTypes>
CameraProjectionModel<DataTypes>::~CameraProjectionModel()
{}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::init()
{
    d_componentState = ComponentState::Valid;
    Inherit1::init();

    if(m_state == nullptr)
    {
        msg_error() << "There is no mechanical state associated with this node.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    internalInit();
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::reinit()
{
    internalInit();
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::internalInit()
{
    if(!d_directions.isSet())
    {
        setDefaultDirections();
    }
    else
    {
        normalizeDirections();
    }

    if(!d_useDirections.isSet())
    {
        setDefaultUseDirections();
    }
    else
    {
        const auto useDirections = sofa::helper::getReadAccessor(d_useDirections);
        if (std::find(useDirections.begin(), useDirections.end(), true) == useDirections.end())
        {
            setDefaultUseDirections();
            msg_warning(this) << "No direction given in useDirections. Set default to all true.";
        }
    }

    if(!d_indices.isSet())
    {
        msg_warning(this) << "Indices not defined. Default value assigned (0).";
        setIndicesDefaultValue();
    }

    if(d_indices.getValue().size() > m_state->getSize())
    {
        msg_warning(this) << "Indices size cannot be larger than the number of points in state. Launching resize process.";
        resizeIndicesRegardingState();
    }

    if(d_indices.getValue().empty())
    {
        msg_error(this) << "Indices size is zero. Component disabled.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    checkIndicesRegardingState();
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::checkIndicesRegardingState()
{
    ReadAccessor<sofa::Data<VecCoord>> positions = m_state->readPositions();
    const auto& indices = d_indices.getValue();

    if (indices.size() > positions.size())
    {
        msg_error(this) << "Indices size (" << indices.size() 
                        << ") is larger than mechanicalState size (" << positions.size() << ").";
        d_componentState = ComponentState::Invalid;
        return;
    }

    for (unsigned int i = 0; i < indices.size(); ++i)
    {
        if (indices[i] >= positions.size())
        {
            msg_error(this) << "Index at entry " << i << " (" << indices[i] 
                            << ") exceeds mechanicalState position size (" << positions.size() << ").";
            d_componentState = ComponentState::Invalid;
            return;
        }
    }
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> defaultIndices = d_indices;
    defaultIndices.resize(1);
    defaultIndices[0] = 0;
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}

// =============================================================================
// COMPUTE PROJECTED ELLIPSE
// =============================================================================
template<class DataTypes>
Eigen::Matrix<double, 5, 1> CameraProjectionModel<DataTypes>::computeProjectedEllipse(const Coord& pose3D)
{
    Eigen::Vector3d center3D(pose3D[0], pose3D[1], pose3D[2]);
    Eigen::Quaterniond orientation(pose3D[6], pose3D[3], pose3D[4], pose3D[5]); // (w, x, y, z)
    Eigen::Matrix3d R = orientation.toRotationMatrix();

    const double fx = d_focalLength.getValue()[0];
    const double fy = d_focalLength.getValue()[1];
    const double u0 = d_principalPoint.getValue()[0];
    const double v0 = d_principalPoint.getValue()[1];
    const double radius = d_ellipseRadius.getValue();

    // 3D local edge points along local X and Y axes
    Eigen::Vector3d edgeUCam = R * Eigen::Vector3d(radius, 0.0, 0.0) + center3D;
    Eigen::Vector3d edgeVCam = R * Eigen::Vector3d(0.0, radius, 0.0) + center3D;

    double Zc = std::max(center3D.z(), 1e-4);
    double Zu = std::max(edgeUCam.z(), 1e-4);
    double Zv = std::max(edgeVCam.z(), 1e-4);

    // Pin-hole perspective projection
    Eigen::Vector2d p_center(fx * (center3D.x() / Zc) + u0, fy * (center3D.y() / Zc) + v0);
    Eigen::Vector2d p_u(fx * (edgeUCam.x() / Zu) + u0,       fy * (edgeUCam.y() / Zu) + v0);
    Eigen::Vector2d p_v(fx * (edgeVCam.x() / Zv) + u0,       fy * (edgeVCam.y() / Zv) + v0);

    Eigen::Vector2d u = p_u - p_center;
    Eigen::Vector2d v = p_v - p_center;

    Eigen::Matrix<double, 5, 1> ellipse;
    ellipse << p_center.x(), p_center.y(), 2.0 * u.norm(), 2.0 * v.norm(), std::atan2(u.y(), u.x()) * (180.0 / M_PI);
    return ellipse;
}

// =============================================================================
// JACOBIAN
// =============================================================================
template<class DataTypes>
void CameraProjectionModel<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
                                                              DataMatrixDeriv &cMatrix,
                                                              unsigned int &cIndex,
                                                              const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    d_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
    MatrixDeriv& column = *cMatrix.beginEdit();
    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    sofa::Index sizeIndices = indices.size();
    
    auto readAccessor = sofa::helper::getReadAccessor(x);

    const double h = 1e-5; // Finite difference step size

    for (unsigned int i = 0; i < sizeIndices; ++i)
    {
        Coord coord = readAccessor[indices[i]];

        // 1. Nominal / Baseline state
        Eigen::Matrix<double, 5, 1> ellipse_0 = computeProjectedEllipse(coord);

        // 2. Derivatives with respect to Translation (x, y, z)
        Coord coord_x = coord; coord_x[0] += h;
        Coord coord_y = coord; coord_y[1] += h;
        Coord coord_z = coord; coord_z[2] += h;

        Eigen::Matrix<double, 5, 1> dE_dx = (computeProjectedEllipse(coord_x) - ellipse_0) / h;
        Eigen::Matrix<double, 5, 1> dE_dy = (computeProjectedEllipse(coord_y) - ellipse_0) / h;
        Eigen::Matrix<double, 5, 1> dE_dz = (computeProjectedEllipse(coord_z) - ellipse_0) / h;

        // 3. Derivatives with respect to Rotation / Quaternions (qx, qy, qz)
        Coord coord_qx = coord; coord_qx[3] += h;
        Coord coord_qy = coord; coord_qy[4] += h;
        Coord coord_qz = coord; coord_qz[5] += h;

        Eigen::Matrix<double, 5, 1> dE_dqx = (computeProjectedEllipse(coord_qx) - ellipse_0) / h;
        Eigen::Matrix<double, 5, 1> dE_dqy = (computeProjectedEllipse(coord_qy) - ellipse_0) / h;
        Eigen::Matrix<double, 5, 1> dE_dqz = (computeProjectedEllipse(coord_qz) - ellipse_0) / h;

        // 4. Construct the 5 Jacobian rows using Deriv (natively compatible with Rigid3d)
        typename DataTypes::Deriv J_cx(dE_dx[0], dE_dy[0], dE_dz[0], dE_dqx[0], dE_dqy[0], dE_dqz[0]);
        typename DataTypes::Deriv J_cy(dE_dx[1], dE_dy[1], dE_dz[1], dE_dqx[1], dE_dqy[1], dE_dqz[1]);
        typename DataTypes::Deriv J_major(dE_dx[2], dE_dy[2], dE_dz[2], dE_dqx[2], dE_dqy[2], dE_dqz[2]);
        typename DataTypes::Deriv J_minor(dE_dx[3], dE_dy[3], dE_dz[3], dE_dqx[3], dE_dqy[3], dE_dqz[3]);
        typename DataTypes::Deriv J_angle(dE_dx[4], dE_dy[4], dE_dz[4], dE_dqx[4], dE_dqy[4], dE_dqz[4]);

        std::array<typename DataTypes::Deriv, 5> localJacobian = { J_cx, J_cy, J_major, J_minor, J_angle };

        // 5. Fill global SOFA constraint matrix
        for (unsigned int j = 0; j < 5; ++j) 
        {
            MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + (i * 5) + j);
            rowIterator.setCol(indices[i], localJacobian[j]);
        }
    }

    unsigned int totalConstraintRows = sizeIndices * 5;
    cIndex += totalConstraintRows;
    cMatrix.endEdit();

    m_nbLines = cIndex - constraintIndex;
}



template<class DataTypes>
void CameraProjectionModel<DataTypes>::storeResults(vector<double>& delta)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    d_delta.setValue(delta);
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::setDefaultDirections()
{
    VecDeriv directions(Deriv::total_size);
    for (sofa::Size i = 0; i < Deriv::total_size; ++i)
    {
        directions[i][i] = 1.0;
    }
    
    d_directions.setValue(directions);
    d_Jacobian.setValue(directions);
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::setDefaultUseDirections()
{
    Vec<Deriv::total_size, bool> useDirections;
    useDirections.assign(true);
    d_useDirections.setValue(useDirections);
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::normalizeDirections()
{
    WriteAccessor<sofa::Data<VecDeriv>> directions = d_directions;
    directions.resize(Deriv::total_size);
    for (unsigned int i = 0; i < Deriv::total_size; ++i)
    {
        directions[i].normalize();
    }
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::draw(const VisualParams* vparams)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    ReadAccessor<sofa::Data<VecCoord>> positions = m_state->readPositions();
    ReadAccessor<sofa::Data<sofa::type::vector<sofa::Index>>> indices = d_indices;

    vector<Coord> points;
    points.reserve(indices.size());
    for (unsigned int i = 0; i < indices.size(); ++i)
    {
        points.push_back(positions[indices[i]]);
    }

    drawPoints(vparams, points, 10.0f, RGBAColor::green());
}

} // namespace softrobots::constraint
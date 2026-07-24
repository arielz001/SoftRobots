
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
#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <SoftRobots/component/constraint/model/CameraProjectionModel.h>
#include <Eigen/Dense>
//#include <Eigen/Geometry>  // Para cuaterniones
//using namespace std;
using namespace Eigen;

namespace softrobots::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;
using sofa::core::ConstVecCoordId ;
using sofa::helper::WriteAccessor ;
using sofa::helper::ReadAccessor ;
using sofa::type::vector ;
using sofa::type::Vec;
using sofa::type::Vec3;
using sofa::type::RGBAColor;




template<class DataTypes>
CameraProjectionModel<DataTypes>::CameraProjectionModel(MechanicalState* object)
    : Inherit1(object)
    , d_indices(initData(&d_indices, "indices",
                                 "If indices size is lower than target size, \n"
                                 "some target will not be considered"))

    , d_weight(initData(&d_weight, sofa::type::vector<Real>(Deriv::total_size, 1.), "weight",
                          "The parameter sets a weight to the minimization."))

    , d_directions(initData(&d_directions,"directions",
                          "The parameter directions allows to specify the directions in \n"
                          "which you want to solve the position."))

    , d_Jacobian(initData(&d_Jacobian,"JacobianBq",
                          "Jacobian relating node motion to B change\n"
                          "."))
    , d_radiusEllipse(initData(&d_radiusEllipse, "radiusEllipse",
                                 "Radius of the 3D feature ellipse/disk"))
                                 
    , d_focalLength(initData(&d_focalLength, "focalLength",
                                 "Focal length of the camera"))

    , d_principalPoint(initData(&d_principalPoint, "principalPoint",
                                 "Principal point of the camera"))

    , d_useDirections(initData(&d_useDirections,"useDirections",
                              "The parameter useDirections allows to select the directions in \n"
                              "which you want to solve the position. If unspecified, the default \n"
                              "values are all true."))

    , d_delta(initData(&d_delta, "delta","Distance to target"))
{
    d_delta.setReadOnly(true);

    this->addUpdateCallback("updateWeight", {&d_weight}, [this](const sofa::core::DataTracker& t)
                            {
                                SOFA_UNUSED(t);
                                auto weight = sofa::helper::getWriteAccessor(d_weight);
                                if (weight.size() != Deriv::total_size)
                                {
                                    msg_info() << "Wrong size for the data field weight, " << weight.size() <<
                                        " instead of " << Deriv::total_size << ". Resizing, with weight[0] as the default value.";
                                    Real w = weight.empty()? 1.: weight[0];
                                    d_weight.setValue(sofa::type::vector<Real>(Deriv::total_size, w));
                                }
                                return sofa::core::objectmodel::ComponentState::Valid;
                            }, {});
}




template<class DataTypes>
CameraProjectionModel<DataTypes>::~CameraProjectionModel()
{
}





template<class DataTypes>
void CameraProjectionModel<DataTypes>::init()
{
    d_componentState = ComponentState::Valid;
    Inherit1::init();

    if(m_state==nullptr)
    {
        msg_error() << "There is no mechanical state associated with this node. "
                        "the object is deactivated. "
                        "To remove this error message fix your scene possibly by "
                        "adding a MechanicalObject." ;
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
    if(!d_radiusEllipse.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto radiusEllipse = sofa::helper::getReadAccessor(d_radiusEllipse);
    }
    if(!d_focalLength.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        // const auto mum = sofa::helper::getReadAccessor(d_mum);
        const auto focalLength = sofa::helper::getReadAccessor(d_focalLength);
    }
    if(!d_principalPoint.isSet())
    {
        // setDefaultDirections();
    }
    else
    {
        const auto principalPoint = sofa::helper::getReadAccessor(d_principalPoint);
    }



// ############################
    if(!d_directions.isSet())
        setDefaultDirections();
    else
        normalizeDirections();
// ###########################


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
            msg_warning(this) << "No direction given in useDirection. Set default all.";
        }
    }
// #########################
    if(!d_indices.isSet())
    {
        msg_warning(this) <<"Indices not defined. Default value assigned 0.";
        setIndicesDefaultValue();
    }

    if(d_indices.getValue().size() > m_state->getSize())
    {
        msg_warning(this) <<"Indices size can not be larger than the number of point in the context. Launch resize process.";
        resizeIndicesRegardingState();
    }

    if(d_indices.getValue().size() == 0)
    {
        msg_error(this) <<"Indices size is zero. The component will not work.";
        d_componentState = ComponentState::Invalid;
        return;
    }

    checkIndicesRegardingState();
}


template<class DataTypes>
void CameraProjectionModel<DataTypes>::checkIndicesRegardingState()
{
    ReadAccessor<sofa::Data<VecCoord> > positions = m_state->readPositions();

    if(d_indices.getValue().size() > positions.size())
    {
        msg_error(this) << "Indices size is larger than mechanicalState size" ;
        d_componentState = ComponentState::Invalid;
        return;
    }

    const auto& indices = d_indices.getValue();
    for(unsigned int i=0; i<indices .size(); i++)
    {
        if (positions.size() <= indices[i])
        {
            msg_error(this) << "Index at index " << i << " is too large regarding mechanicalState [position] size" ;
            d_componentState = ComponentState::Invalid;
            return;
        }
    }
}




template<class DataTypes>
void CameraProjectionModel<DataTypes>::setIndicesDefaultValue()
{
    WriteAccessor<sofa::Data<vector<unsigned int> > > defaultIndices = d_indices;
    defaultIndices.resize(1);
}

template<class DataTypes>
void CameraProjectionModel<DataTypes>::resizeIndicesRegardingState()
{
    WriteAccessor<sofa::Data<vector<unsigned int>>> indices = d_indices;
    indices.resize(m_state->getSize());
}





Eigen::Matrix<double, 5, 1> calculateProjectedEllipse(
    double x, double y, double z,
    const Eigen::Matrix3d& R,
    double radius,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint)
{
    // 1. Centro de la elipse en píxeles (Proyección perspectiva Pinhole)
    double u = focalLength[0] * (x / z) + principalPoint[0];
    double v = focalLength[1] * (y / z) + principalPoint[1];

    // 2. Extraer la normal del disco 3D (tercera columna de la matriz de rotación)
    Eigen::Vector3d normal = R.col(2); 

    // inclination
    double cos_tilt = std::abs(normal(2));
    if (cos_tilt < 1e-3) cos_tilt = 1e-3; // Evitar división por cero si está de canto

    // Semiejes
    double semi_a = focalLength[0] * (radius / z); 
    double semi_b = semi_a * cos_tilt;           

    // 5. Angle
    double alpha = std::atan2(normal(1), normal(0));

    Eigen::Matrix<double, 5, 1> ellipse;
    ellipse << u, v, semi_a, semi_b, alpha;
    return ellipse;
}


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
    
    // Accesores de SOFA para la cámara y el disco
    const auto focalLength = d_focalLength.getValue();
    const auto principalPoint = d_principalPoint.getValue();
    const double radius = d_radiusEllipse.getValue();

    auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
    const auto& readAccessor = sofa::helper::getReadAccessor(x);

    // Epsilon para derivadas numéricas (un valor más fino para coordenadas de pantalla)
    const double Cambio = 1e-4;

    for (const auto& coord : readAccessor) {
        // 1. Obtener Posición 3D
        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // 2. Obtener Orientación 3D (Cuaternión -> Matriz de Rotación)
        Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]); // (w, x, y, z)
        Eigen::Matrix3d R = q.toRotationMatrix();

        // 3. Estado Base de la Elipse 2D [u, v, a, b, alpha]
        Eigen::Matrix<double, 5, 1> E_0 = calculateProjectedEllipse(x_pos, y_pos, z_pos, R, radius, focalLength, principalPoint);

        // --- DERIVADAS RESPECTO A X ---
        Eigen::Matrix<double, 5, 1> E_x = calculateProjectedEllipse(x_pos + Cambio, y_pos, z_pos, R, radius, focalLength, principalPoint);
        Eigen::Matrix<double, 5, 1> dE_dx = (E_x - E_0) / Cambio;

        // --- DERIVADAS RESPECTO A Y ---
        Eigen::Matrix<double, 5, 1> E_y = calculateProjectedEllipse(x_pos, y_pos + Cambio, z_pos, R, radius, focalLength, principalPoint);
        Eigen::Matrix<double, 5, 1> dE_dy = (E_y - E_0) / Cambio;

        // --- DERIVADAS RESPECTO A Z ---
        Eigen::Matrix<double, 5, 1> E_z = calculateProjectedEllipse(x_pos, y_pos, z_pos + Cambio, R, radius, focalLength, principalPoint);
        Eigen::Matrix<double, 5, 1> dE_dz = (E_z - E_0) / Cambio;

        // --- DERIVADAS RESPECTO A ROTACIONES (Perturbación en ejes X, Y, Z de rotación) ---
        // Perturbación Roll (alrededor de X)
        Eigen::Matrix3d R_rx = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitX());
        Eigen::Matrix<double, 5, 1> dE_dRx = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_rx, radius, focalLength, principalPoint) - E_0) / Cambio;

        // Perturbación Pitch (alrededor de Y)
        Eigen::Matrix3d R_ry = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitY());
        Eigen::Matrix<double, 5, 1> dE_dRy = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_ry, radius, focalLength, principalPoint) - E_0) / Cambio;

        // Perturbación Yaw (alrededor de Z)
        Eigen::Matrix3d R_rz = R * Eigen::AngleAxisd(Cambio, Eigen::Vector3d::UnitZ());
        Eigen::Matrix<double, 5, 1> dE_dRz = (calculateProjectedEllipse(x_pos, y_pos, z_pos, R_rz, radius, focalLength, principalPoint) - E_0) / Cambio;

        // 4. Llenar la Matriz Jacobiana (5 Filas: u, v, a, b, alpha)
        for (int i = 0; i < 5; ++i) {
            Jacobian[i] = sofa::type::Vec<6, double>(
                dE_dx[i],   // d/dx
                dE_dy[i],   // d/dy
                dE_dz[i],   // d/dz
                dE_dRx[i],  // d/dRoll
                dE_dRy[i],  // d/dPitch
                dE_dRz[i]   // d/dYaw
            );
        }
    }

    // 5. Escribir las 5 restricciones en el sistema global de SOFA
    unsigned int index = 0;
    for (unsigned j = 0; j < 5; j++) { // 5 filas para los 5 parámetros de la elipse
        MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
        rowIterator.setCol(0, Jacobian[j]);
        index++;
    }

    cIndex += index;
    cMatrix.endEdit();
    m_nbLines = cIndex - constraintIndex;
}





// template<class DataTypes>
// void CameraProjectionModel<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
//                                                         DataMatrixDeriv &cMatrix,
//                                                         unsigned int &cIndex,
//                                                         const DataVecCoord &x)
// {
//     if(d_componentState.getValue() != ComponentState::Valid)
//         return;

//     SOFA_UNUSED(cParams);
//     SOFA_UNUSED(x);





//     d_constraintIndex.setValue(cIndex);
//     const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);
//     MatrixDeriv& column = *cMatrix.beginEdit();
//     const auto& indices = sofa::helper::getReadAccessor(d_indices);
//     sofa::Index sizeIndices = indices.size();
//     const auto& useDirections = sofa::helper::getReadAccessor(d_useDirections);
//     auto directions = sofa::helper::getWriteAccessor(d_directions);
//     auto Jacobian = sofa::helper::getWriteAccessor(d_Jacobian);
//     auto PosSensor = sofa::helper::getWriteAccessor(d_PosSensor);
//     auto mum = sofa::helper::getWriteAccessor(d_mum);
//     auto weight = sofa::helper::getReadAccessor(d_weight);

// //    std::cout << "datos d epos AAAAAAA: " << x << std::endl;
//     auto& data = x;  // Usamos directamente el ReadAccessor    double dBx_drz = 0;
//     Eigen::Vector3d B_calculada;  // Variable global o de ámbito extendido
//     auto readAccessor = sofa::helper::getReadAccessor(data);  // Obtén un ReadAccessor
//     Eigen::Vector3d B_0;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_1;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_2;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_3;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_Phi;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_Phi_devz;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdR_z;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdR_x;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdR_y;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdTheta;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdPhi;  // Variable global o de ámbito extendido
//     Eigen::Vector3d dBdPhi_devz;  // Variable global o de ámbito extendido

// //    Eigen::Vector3d dBdR_y;  // Variable global o de ámbito extendido
//     const auto Cambio = 0.01;
//     // PosSensor[0] =  sofa::type::Vec<6, double>(0,0,0,0,0,24); //PosSensor se define aqui, debo linkearlo con la escena o algo así, 0,1,2 Para sensor, 3,4,5 Para Iman
//     // PosSensor[0] =  sofa::type::Vec<6, double>(2,3,4,5,6,7); //PosSensor se define aqui, debo linkearlo con la escena o algo así
//     // PosSensor[1] =  sofa::type::Vec<6, double>(8,9,10,11,12,13); //PosSensor se define aqui, debo linkearlo con la escena o algo así


//     for (const auto& coord : readAccessor) {
//             // Acceder a todas las componentes de cada Vec<2, double> y mostrar las tres componentes si es posible
//             // msg_warning() << "Coord: (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";  // Asumiendo que 'Vec<2, double>' tiene tres componentes
//             Eigen::Quaterniond MiR(coord[6], coord[3], coord[4], coord[5]);  // (w, x, y, z)
//             // std::cout << "MiR de CameraProjectionModel.inl : " << MiR.coeffs().transpose() << "\n";
//             Eigen::Matrix3d rotation_matrix = MiR.toRotationMatrix();
//             // std::cout << "Matriz de rotación:\n" << rotation_matrix << std::endl;
//             // Extraer la última columna
//             const double mu_x = rotation_matrix(0, 2);  // Elemento (0, 2)
//             const double mu_y = rotation_matrix(1, 2);  // Elemento (1, 2)
//             const double mu_z = rotation_matrix(2, 2);  // Elemento (2, 2)
//             // std::cout << "mu_x: " << mu_x << ", mu_y: " << mu_y << ", mu_z: " << mu_z << std::endl;
//             const double ajuste_x = PosSensor[0][0] ;
//             const double ajuste_y = PosSensor[0][1];
//             const double ajuste_z = PosSensor[0][2];
//             // std::cout << "PosSensor[0] de CameraProjectionModel.inl : " << PosSensor[0] << "\n";
//     //        const double ajuste_z = 2.7 - 15; // Este ajuste era para corroborar calculo en c++ con el de python
//             B_calculada = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2],mum[0][0],mu_x,mu_y,mu_z);
//         //    std::cout << "Campo magnético B_c alculado c++ AAAAAAAAAAAAAA: " << B_calculada.transpose() << std::endl;
//             B_0 = Calculo_B_Test(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2],mum[0][0],mu_x,mu_y,mu_z);
//             // dB/drX
//             B_2 = Calculo_B_Test(coord[0] - ajuste_x + Cambio,coord[1]- ajuste_y,coord[2],mum[0][0],mu_x,mu_y,mu_z);
//             dBdR_x = (B_2 - B_0)/Cambio;
//             // std::cout << "dBdR_x" << dBdR_x.transpose() << std::endl;
//             // dB/drY
//             B_3 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y+ Cambio,coord[2],mum[0][0],mu_x,mu_y,mu_z);
//             dBdR_y = (B_3 - B_0)/Cambio;
//             // std::cout << "dBdR_y" << dBdR_y.transpose() << std::endl;
//             // dB/drZ
//             B_1 = Calculo_B_Test(coord[0] - ajuste_x,coord[1]- ajuste_y,coord[2]    + Cambio    ,mum[0][0],mu_x,mu_y,mu_z);
//             dBdR_z = (B_1 - B_0)/Cambio;
//             // std::cout << "dBdR_z" << dBdR_z.transpose() << std::endl;

            
// //            ---------------- Angulos :-----------------------------------------------------------

// //            Eigen::Vector3d Mu(0.44229157 ,-0.32357449, -0.8364674); // Ejemplo de entrada, Resultado: ThetaRecovered (rad): 2.81208787) PhiRecovered: (rad) -0.4863910100000001
//             Eigen::Vector3d Mu(mu_x, mu_y, mu_z );
//             auto [theta, phi] = recoverThetaAndPhi(Mu);
//             // std::cout << "Theta: " << theta << "\n";
//             // std::cout << "Phi: " << phi << "\n";

//             // Crear rotaciones alrededor de X y Y
//             Eigen::AngleAxisd Rx(theta + Cambio, Eigen::Vector3d::UnitX());
//             Eigen::AngleAxisd Ry(phi, Eigen::Vector3d::UnitY());
//             Eigen::Quaterniond MiR_2 = Ry * Rx;
//             Eigen::Vector3d Mu_theta = MiR_2 * Eigen::Vector3d(0, 0, 1);
//             // std::cout << "Mu_theta: " << Mu_theta.transpose() << std::endl;

// //            auto [theta_2, phi_2] = recoverThetaAndPhi(Mu_theta);
// //            std::cout << "Theta_2: " << theta_2 << "\n";
// //            std::cout << "Phi_2: " << phi_2 << "\n";


//             // db/dTheta
//             B_1 = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_theta[0],Mu_theta[1],Mu_theta[2]);
//             dBdTheta= (B_1 - B_0)/ (Cambio);
//             // std::cout << "dBdTheta : " << dBdTheta.transpose() << std::endl;

//             // db/dPhi
//             Eigen::AngleAxisd Rx_Phi(theta, Eigen::Vector3d::UnitX());
//             Eigen::AngleAxisd Ry_Phi(phi+ Cambio, Eigen::Vector3d::UnitY());
//             Eigen::Quaterniond MiR_Phi = Ry_Phi * Rx_Phi;
//             Eigen::Vector3d Mu_Phi = MiR_Phi* Eigen::Vector3d(0, 0, 1);
//             B_Phi = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_Phi[0],Mu_Phi[1],Mu_Phi[2]);
//             dBdPhi= (B_Phi - B_0)/ (Cambio);
//             // std::cout << "dBdPhi : " << dBdPhi.transpose() << std::endl;

// //            ---------------- Angulo para eje z :-----------------------------------------------------------
//             // Eigen::Vector3d Mu(mu_x, mu_y, mu_z );
//             auto [theta_devz, phi_devz] = recoverThetaAndPhi_DerivadaZ(Mu);
//             // std::cout << "Theta_devz: " << theta_devz << "\n";
//             // std::cout << "Phi_devz: " << phi_devz << "\n";     

//             // Crear rotaciones alrededor de X y Y
//             Eigen::AngleAxisd Rx_devz(theta_devz , Eigen::Vector3d::UnitX());
//             Eigen::AngleAxisd Rz_devz(phi_devz + Cambio, Eigen::Vector3d::UnitZ());
//             Eigen::Quaterniond MiR_devz = Rz_devz * Rx_devz;
//             Eigen::Vector3d Mu_devz = MiR_devz * Eigen::Vector3d(0, 0, 1);

//             B_Phi_devz = Calculo_B_Test(coord[0] - ajuste_x ,coord[1]- ajuste_y, coord[2],mum[0][0], Mu_devz[0],Mu_devz[1],Mu_devz[2]);
//             // std::cout << "mum de CameraProjectionModel.inl : " << mum[0][0] << "\n";     
//             dBdPhi_devz= (B_Phi_devz - B_0)/ (Cambio);
//         }


//     double dBx_drx = dBdR_x[0];
//     double dBy_drx = dBdR_x[1];
//     double dBz_drx = dBdR_x[2];
//     double dBx_dry = dBdR_y[0];
//     double dBy_dry = dBdR_y[1];
//     double dBz_dry = dBdR_y[2];
//     double dBx_drz = dBdR_z[0];
//     double dBy_drz = dBdR_z[1];
//     double dBz_drz = dBdR_z[2];

//     double dBx_dTheta = dBdTheta[0];
//     double dBy_dTheta = dBdTheta[1];
//     double dBz_dTheta = dBdTheta[2];

//     double dBx_dPhi = dBdPhi[0];
//     double dBy_dPhi = dBdPhi[1];
//     double dBz_dPhi = dBdPhi[2];

//     double dBx_dPhi_devz = dBdPhi_devz[0];
//     double dBy_dPhi_devz = dBdPhi_devz[1];
//     double dBz_dPhi_devz = dBdPhi_devz[2];

// //    Jacobian.push_back(VecDeriv(0.0, 0.0, dBx_drz, 0.0, 0.0, 0.0));


// //    Jacobian[0] = sofa::type::Vec<6, double>(17114, 0.0, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
// //    Jacobian[1] = sofa::type::Vec<6, double>(0.0, 17114, 0.0, 0.0, 0.0, 0.0);  // Asignar un nuevo valor
// //    Jacobian[2] = sofa::type::Vec<6, double>(-126.7, -126.7, -33978, 0.0, 0.0, 0.0);  // Asignar un nuevo valor

//         Jacobian[0] = sofa::type::Vec<6, double>(dBx_drx, dBx_dry, dBx_drz,     dBx_dTheta, dBx_dPhi, dBx_dPhi_devz);  // Asignar un nuevo valor
//         Jacobian[1] = sofa::type::Vec<6, double>(dBy_drx, dBy_dry, dBy_drz,     dBy_dTheta, dBy_dPhi, dBy_dPhi_devz);  // Asignar un nuevo valor
//         Jacobian[2] = sofa::type::Vec<6, double>(dBz_drx, dBz_dry, dBz_drz,     dBz_dTheta, dBz_dPhi, dBz_dPhi_devz);  // Asignar un nuevo valor



//     unsigned int index = 0;

//     for (unsigned j = 0; j < 3; j++) {
//         MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex + index);
//     //    std::cout << "datos de directions: " << directions[j] << std::endl;

//         rowIterator.setCol(0, Jacobian[j]);
//         index++;
//     }


//     // for (unsigned i=0; i<sizeIndices; i++)
//     // {
//     //     for(sofa::Size j=0; j<Deriv::total_size; j++)
//     //     {
//     //         if(useDirections[j])
//     //         {
//     //             MatrixDerivRowIterator rowIterator = column.writeLine(constraintIndex+index);
//     //             std::cout << "Tipo de directions: " << typeid(directions).name() << std::endl;
//     //             rowIterator.setCol(indices[i], directions[j]);
//     //             index++;
//     //         }
//     //     }
//     // }

//     cIndex += index;
//     cMatrix.endEdit();
    
//     m_nbLines = cIndex - constraintIndex;
// }


template<class DataTypes>


void CameraProjectionModel<DataTypes>::storeResults(vector<double> &delta)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    d_delta.setValue(delta);
}


template<class DataTypes>
void CameraProjectionModel<DataTypes>::setDefaultDirections()
{
    VecDeriv directions(Deriv::total_size);
    for(sofa::Size i=0; i<Deriv::total_size; i++)
        directions[i][i] = 1.;
    d_directions.setValue(directions);
    d_Jacobian.setValue(directions);
    // d_PosSensor.setValue(directions);
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
    for(unsigned int i=0; i<Deriv::total_size; i++)
        directions[i].normalize();
}


template<class DataTypes>
void CameraProjectionModel<DataTypes>::draw(const VisualParams* vparams)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    ReadAccessor<sofa::Data<VecCoord> > positions = m_state->readPositions();
    ReadAccessor<sofa::Data<sofa::type::vector<sofa::Index>> > indices = d_indices;
    vector<Coord> points;
    points.reserve(indices.size());
    for (unsigned int i=0; i<indices.size(); i++)
    {
        points.push_back(positions[indices[i]]);
    }
    drawPoints(vparams, points, 10.0f, RGBAColor::green());
}

} // namespace

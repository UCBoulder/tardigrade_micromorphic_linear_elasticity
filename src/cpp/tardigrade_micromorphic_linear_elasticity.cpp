/*
 * tardigrade_micromorphic_linear_elasticity.cpp
 *
 * An implimentation of linear elasticity in the micromorphic context.
 *
 * Based around a quadratic form of the Helmholtz free energy:
 * \f$\rho \psi = \frac{1}{2} E_{IJ} A_{IJKL} E_{KL} + \frac{1}{2} \mathcal{E}_{IJ} B_{IJKL} \mathcal{E}_{KL} 
 *              + \frac{1}{2} \Gamma_{IJK} C_{IJKLMN} \Gamma_{LMN} + E_{IJ} D_{IJKL} \mathcal{E}_{KL}\f$
 */

#include<tardigrade_hydraMicromorphicLinearElasticity.h>
#include<tardigrade_micromorphic_linear_elasticity.h>

namespace tardigradeMicromorphicLinearElasticity{

    void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &PK2Stress: The second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         */

        //Compute the required deformation measures
        variableVector RCG, Psi, Gamma;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                  RCG, Psi, Gamma ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                                                PK2Stress, referenceMicroStress,
                                                                                referenceHigherOrderStress ) );

        return;
    }

    void linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation,
                                                       const variableVector &Psi, const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
         * \param &Psi: The micro-deformation measure Psi
         * \param &Gamma: The higher order deformation measure Gamma
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &PK2Stress: The second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        variableVector invRCG = tardigradeVectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the strain measures
        variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
        variableVector microStrain   = Psi - eye;

        //Compute the higher order stress
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress ) );

        //Compute the first common term for the PK2 and symmetric micro-stress
        variableVector term1;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 ) );

        //Compute the second common term for the PK2 and symmetric micro-stress
        variableVector invRCGPsi;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi ) );

        variableVector term2;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2 ) );

        //Compute the third common term for the PK2 and symmetric micro-stress
        variableVector invRCGGamma;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma ) );

        variableVector term3;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3 ) );

        //Construct the PK2 and reference symmetric stresses
        PK2Stress            = term1 + term2 + term3;

        variableVector symmTerm2Term3;
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3 ) );
        referenceMicroStress = term1 + 2 * symmTerm2Term3;

        return;
    }

    void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdChi, variableMatrix &dPK2StressdGradChi,
                                        variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdChi,
                                        variableMatrix &dReferenceMicroStressdGradChi, variableMatrix &dMdF, variableMatrix &dMdGradChi ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * Also computes the Jacobians
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &PK2Stress: The second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         * \param &dPK2StressdF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
         * \param &dPK2StressdChi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
         * \param &dPK2StressdGradChi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
         * \param &dReferenceMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient.
         * \param &dReferenceMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
         * \param &dReferenceMicroStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the micro deformation.
         * \param &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
         * \param &dMdGradChi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
         */

        //Assume 3d
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        //Compute the required deformation measures
        variableVector RCG, Psi, Gamma;
        variableMatrix dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                     RCG, Psi, Gamma, dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi ) );

        variableMatrix dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma;
        variableMatrix dReferenceMicroStressdRCG, dReferenceMicroStressdPsi, dReferenceMicroStressdGamma;
        variableMatrix dMdGamma;

        TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                          PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                          dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma,
                                                          dReferenceMicroStressdRCG, dReferenceMicroStressdPsi,
                                                          dReferenceMicroStressdGamma, dMdGamma ) );

        dPK2StressdF = tardigradeVectorTools::dot( dPK2StressdRCG, dRCGdF )
                     + tardigradeVectorTools::dot( dPK2StressdPsi, dPsidF )
                     + tardigradeVectorTools::dot( dPK2StressdGamma, dGammadF );

        dPK2StressdChi = tardigradeVectorTools::dot( dPK2StressdPsi, dPsidChi );

        dPK2StressdGradChi = tardigradeVectorTools::dot( dPK2StressdGamma, dGammadGradChi );

        dReferenceMicroStressdF = tardigradeVectorTools::dot( dReferenceMicroStressdRCG, dRCGdF )
                                + tardigradeVectorTools::dot( dReferenceMicroStressdPsi, dPsidF )
                                + tardigradeVectorTools::dot( dReferenceMicroStressdGamma, dGammadF );

        dReferenceMicroStressdChi = tardigradeVectorTools::dot( dReferenceMicroStressdPsi, dPsidChi );

        dReferenceMicroStressdGradChi = tardigradeVectorTools::dot( dReferenceMicroStressdGamma, dGammadGradChi );

        dMdF = tardigradeVectorTools::dot( dMdGamma, dGammadF );
        dMdGradChi = tardigradeVectorTools::dot( dMdGamma, dGammadGradChi );

        return;
    }

    void linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                       const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress,
                                                       variableMatrix &dPK2StressdRCG, variableMatrix &dPK2StressdPsi,
                                                       variableMatrix &dPK2StressdGamma,
                                                       variableMatrix &dReferenceMicroStressdRCG,
                                                       variableMatrix &dReferenceMicroStressdPsi,
                                                       variableMatrix &dReferenceMicroStressdGamma,
                                                       variableMatrix &dMdGamma ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
         * \param &Psi: The micro-deformation measure Psi
         * \param &Gamma: The higher order deformation measure Gamma
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &PK2Stress: The second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         * \param &dPK2StressdRCG: The Jacobian of the PK2 stress w.r.t. the right Cauchy-Green
         *     deformation metric.
         * \param &dPK2StressdPsi: The Jacobian of the PK2 stress w.r.t. the micro deformation 
         *     metric.
         * \param &dPK2StressdGamma: The Jacobian of the PK2 stress w.r.t. the higher order 
         *     deformation measure.
         * \param &dReferenceMicroStressdRCG: The Jacobian of the reference micro stress w.r.t. the 
         *     right Cacuhy-Green deformation metric.
         * \param &dReferenceMicroStressdPsi: The Jacobian of the reference micro stress w.r.t. the 
         *     micro deformation measure.
         * \param &dReferenceMicroStressdGamma: The Jacobian of the reference micro stress w.r.t. the 
         *     higher order deformation measure.
         * \param &dMdGamma: The Jacobian of the reference higher order stress w.r.t. 
         *     the higher order deformation measure.
         */

        //Assume 3d
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        variableVector invRCG = tardigradeVectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the strain measures
        variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
        variableVector microStrain   = Psi - eye;

        //Compute the higher order stress
        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress, dMdGamma ) );

        //Compute the first common term for the PK2 and symmetric micro-stress
        variableVector term1;

        variableMatrix dTerm1dRCG, dTerm1dPsi;
        TARDIGRADE_ERROR_TOOLS_CATCH(computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1,
                                           dTerm1dRCG, dTerm1dPsi ) );

        //Assemble term1 jacobians w.r.t. F and Chi
        dTerm1dRCG *= 0.5;

        //Compute the second common term for the PK2 and symmetric micro-stress
        variableVector invRCGPsi;
        variableMatrix dInvRCGPsidRCG, dInvRCGPsidPsi;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi ) );

        variableVector term2;
        variableMatrix dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi;
        TARDIGRADE_ERROR_TOOLS_CATCH(computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                                                dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi ) );

        dTerm2dRCG *= 0.5;
        dTerm2dRCG += tardigradeVectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidRCG );

        dTerm2dPsi += tardigradeVectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidPsi );

        //Compute the third common term for the PK2 and symmetric micro-stress
        variableVector invRCGGamma;
        variableMatrix dInvRCGGammadRCG, dInvRCGGammadGamma;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma ) );

        variableVector term3;
        variableMatrix dTerm3dInvRCGGamma, dTerm3dM;
        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvRCGGamma, dTerm3dM ) );

        variableMatrix dTerm3dRCG = tardigradeVectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadRCG );
        variableMatrix dTerm3dGamma = tardigradeVectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadGamma )
                                    + tardigradeVectorTools::dot( dTerm3dM, dMdGamma );

        //Construct the PK2 and reference symmetric stresses
        PK2Stress            = term1 + term2 + term3;

        dPK2StressdRCG    = dTerm1dRCG + dTerm2dRCG + dTerm3dRCG;
        dPK2StressdPsi    = dTerm1dPsi + dTerm2dPsi;
        dPK2StressdGamma  = dTerm3dGamma;

        variableVector symmTerm2Term3;
        variableMatrix dSymmTerm2Term3dTerm2Term3;
        tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 );
        referenceMicroStress = term1 + 2 * symmTerm2Term3;

        dReferenceMicroStressdRCG = dTerm1dRCG + 2 * ( tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dRCG )
                                                     + tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dRCG ) );

        dReferenceMicroStressdPsi = dTerm1dPsi + 2 * tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dPsi );
        dReferenceMicroStressdGamma = 2 * tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dGamma );

        return;
    }

    void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma ){
        /*!
         * Compute the deformation measures
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         * \f$\Psi_{IJ} = F_{iI} \Chi_{iJ}\f$
         * \f$\Gamma_{IJK} = F_{iI} \Chi_{iJ, K}\f$
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &rightCauchyGreen: The Right Cauchy-Green deformation tensor
         * \param &Psi: The micro-deformation measure
         * \param &Gamma: The gradient micro-deformation measure
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma ) );

        return;

    }

    void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                         variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidChi,
                                         variableMatrix &dGammadF, variableMatrix &dGammadGradChi ){
        /*!
         * Compute the deformation measures
         * \f$C_{IJ} = F_{iI} F_{iJ}\f$
         * \f$\Psi_{IJ} = F_{iI} \Chi_{iJ}\f$
         * \f$\Gamma_{IJK} = F_{iI} \Chi_{iJ, K}\f$
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &rightCauchyGreen: The Right Cauchy-Green deformation tensor
         * \param &Psi: The micro-deformation measure
         * \param &Gamma: The gradient micro-deformation measure
         * \param &dCdF: The gradient of the right Cauchy green deformation tensor w.r.t. 
         *     the deformation gradient.
         * \param &dPsidF: The gradient of Psi w.r.t. the deformation gradient.
         * \param &dPsidChi: The gradient of Psi w.r.t. the microDeformation.
         * \param &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * \param &dGammadGradChi: The gradient of Gamma w.r.t. the spatial gradient of Chi
         */

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen, dCdF ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi, dPsidF, dPsidChi ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma, dGammadF, dGammadGradChi ) );

        return;
    }

    void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain, 
                                        const parameterVector &A, const parameterVector &D, variableVector &term1 ){
        /*!
         * Compute the first term for the linear elastic model
         * \f$term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL}\f$
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &microStrain: The micro-strain
         * \param &A: The A stiffness matrix
         * \param &D: The D stiffness matrix
         * \param &term1: The first term.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( greenLagrangeStrain.size() == dim * dim, "The green lagrange strain must have a length of 9" );

        TARDIGRADE_ERROR_TOOLS_CHECK( microStrain.size() == dim * dim, "The micro-strain must have a length of 9" );

        TARDIGRADE_ERROR_TOOLS_CHECK( A.size() == dim * dim * dim * dim, "A must have a size of 3**4" );

        TARDIGRADE_ERROR_TOOLS_CHECK( D.size() == dim * dim * dim * dim, "D must have a size of 3**4" );

        //Compute the first common term for the PK2 and symmetric micro-stress
        term1 = variableVector( dim * dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        term1[ dim * I + J ] += A[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * greenLagrangeStrain[ dim * K + L ]
                                              + D[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * microStrain[ dim * K + L ];
                    }
                }
            }
        }

        return;
    }

    void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain, 
                                        const parameterVector &A, const parameterVector &D, variableVector &term1,
                                        variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain ){
        /*!
         * Compute the first term for the linear elastic model
         * \f$term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL}\f$
         *
         * Also return the Jacobian
         * \f$\frac{\partial term^1_{IJ} }{ E_{MN} } = A_{IJMN}\f$
         * \f$\frac{\partial term^1_{IJ} }{ \mathcal{E}_{MN} } = D_{IJMN}\f$
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain.
         * \param &microStrain: The micro-strain
         * \param &A: The A stiffness matrix
         * \param &D: The D stiffness matrix
         * \param &term1: The first term.
         * \param &dTerm1dGreenLagrangeStrain: The derivative of term1 w.r.t. the Green-Lagrange strain
         * \param &dTerm1dMicroStrain: The derivative of term1 w.r.t. the micro strain
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 ) );

        //Compute the first common term for the PK2 and symmetric micro-stress
        dTerm1dGreenLagrangeStrain = variableMatrix( term1.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
        dTerm1dMicroStrain = variableMatrix( term1.size(), variableVector( microStrain.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        dTerm1dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] = A[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                        dTerm1dMicroStrain[ dim * I + J ][ dim * M + N ] = D[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                    }
                }
            }
        }

        return;
    }

    void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2 ){
        /*!
         * Compute the second term from the linear elastic constitutive model
         *
         * \f$term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}\f$
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
         * \param &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
         * \param &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
         * \param &B: The B stiffness matrix
         * \param &D: The D stiffness matrix
         * \param &term2: The second term.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( greenLagrangeStrain.size() == dim * dim, "The green lagrange strain must have a length of 9" );

        TARDIGRADE_ERROR_TOOLS_CHECK( microStrain.size() == dim * dim, "The micro-strain must have a length of 9" );

        TARDIGRADE_ERROR_TOOLS_CHECK( invCPsi.size() == dim * dim, "invCPsi must have a size of 9" );

        TARDIGRADE_ERROR_TOOLS_CHECK( B.size() == dim * dim * dim * dim, "B must have a size of 3**4" );

        TARDIGRADE_ERROR_TOOLS_CHECK( D.size() == dim * dim * dim * dim, "D must have a size of 3**4" );

        term2 = variableVector( greenLagrangeStrain.size(), 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int Q = 0; Q < dim; Q++ ){
                            term2[ dim * I + J] += ( B[ dim * dim * dim * I + dim * dim * Q + dim * K + L ] * microStrain[ dim * K + L ]
                                                 + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + Q ] )
                                                 * invCPsi[ dim * J + Q ];
                        }
                    }
                }
            }
        }

        return;
    }

    void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                        variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi ){
        /*!
         * Compute the second term from the linear elastic constitutive model
         *
         * \f$term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}\f$
         *
         * Also return the Jacobians
         * \f$\frac{ \partial term^2_{IJ} }{ \partial E_{MN} } = D_{MNIK} C_{JR}^{-1} \Psi_{RK}\f$
         * \f$\frac{ \partial term^2_{IJ} }{ \partial \mathcal{E}_{MN} } = B_{IKMN} C_{JR}^{-1} \Psi_{RK}\f$
         * \f$\frac{ \partial term^2_{IJ} }{ \partial C_{MO}^{-1} \Psi_{ON} } = \left( B_{INKL} \mathcal{E}_{KL} + E_{KL} D_{KLIN} \right) \delta_{JM}\f$
         *
         * \param &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
         * \param &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
         * \param &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
         * \param &B: The B stiffness matrix
         * \param &D: The D stiffness matrix
         * \param &term2: The second term.
         * \param &dTerm2dGreenLagrangeStrain: The jacobian of term 2 w.r.t. the Green-Lagrange strain.
         * \param &dTerm2dMicroStrain: The jacobian of term 2 w.r.t. the microStrain.
         * \param &dTerm2dInvCPsi: The jacobian of term 2 w.r.t. \f$C_{JR}^{-1} \Psi_{RQ}\f$
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invCPsi, B, D, term2 ) );

        //Compute the Jacobians
        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );
        dTerm2dGreenLagrangeStrain = variableMatrix( term2.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
        dTerm2dMicroStrain         = variableMatrix( term2.size(), variableVector( microStrain.size(), 0 ) );
        dTerm2dInvCPsi             = variableMatrix( term2.size(), variableVector( invCPsi.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int K = 0; K < dim; K++ ){
                            dTerm2dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] += D[ dim * dim * dim * M + dim * dim * N + dim * I + K] * invCPsi[ dim * J + K ];
                            dTerm2dMicroStrain[ dim * I + J ][ dim * M + N ] += B[ dim * dim * dim * I + dim * dim * K + dim * M + N] * invCPsi[ dim * J + K ];
                            for ( unsigned int L = 0; L < dim; L++ ){
                                dTerm2dInvCPsi[ dim * I + J ][ dim * M + N ] += ( B[ dim * dim * dim * I + dim * dim * N + dim * K + L ] * microStrain[ dim * K + L ] + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + N ] ) * eye[ dim * J + M ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the higher order stress in the reference configuration.
         * \f$M_{IJK} = C_{JKILMN} \Gamma_{LMN}\f$
         *
         * \param &Gamma: The micro-gradient deformation measure.
         * \param &C: The C stiffness tensor.
         * \param &referenceHigherOrderStress: The higher order stress in the reference 
         *     configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( Gamma.size() == dim * dim * dim, "Gamma must have a length of 27" );

        TARDIGRADE_ERROR_TOOLS_CHECK( C.size() == dim * dim * dim * dim * dim * dim, "The C stiffness tensor have a length of 3**6" );

        referenceHigherOrderStress = variableVector( dim * dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                referenceHigherOrderStress[ dim * dim * I + dim * J + K ] += C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * L + dim * M + N ] * Gamma[ dim * dim * L + dim * M + N ];
                            }
                        }
                    }
                }
            }
        }
        return;
    }

    void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress,
                                                variableMatrix &dReferenceHigherOrderStressdGamma ){
        /*!
         * Compute the higher order stress in the reference configuration.
         * \f$M_{IJK} = C_{JKILMN} Gamma_{LMN}\f$
         *
         * Also compute the Jacobian
         * \f$\frac{ \partial M_{IJK} }{\partial \Gamma_{OPQ} } = C_{JKIOPQ}\f$
         *
         * \param &Gamma: The micro-gradient deformation measure.
         * \param &C: The C stiffness tensor.
         * \param &referenceHigherOrderStress: The higher order stress in the reference 
         *     configuration.
         * \param &dReferenceHigherOrderStressdGamma: The derivative of the higher order stress
         *     in the reference configuration w.r.t. Gamma
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress ) );

        //Assemble the Jacobian
        dReferenceHigherOrderStressdGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int O = 0; O < dim; O++ ){
                        for ( unsigned int P = 0; P < dim; P++ ){
                            for ( unsigned int Q = 0; Q < dim; Q++ ){
                                dReferenceHigherOrderStressdGamma[ dim * dim * I + dim * J + K ][ dim * dim * O + dim * P + Q ] +=
                                    C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * O + dim * P + Q ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void computeLinearElasticTerm3( const variableVector &invCGamma, 
                                        const variableVector &referenceHigherOrderStress, variableVector &term3 ){
        /*!
         * Compute the value of the third term in the micromorphic linear elasticity formulation.
         * \f$term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR}\f$
         *
         * \param &invCGamma: \f$C_{JS}^{-1} \Gamma_{SQR}\f$
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &term3: The third term in the linear elastic equation.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( invCGamma.size() == dim * dim * dim, "invCGamma must have a size of 27" );

        TARDIGRADE_ERROR_TOOLS_CHECK( referenceHigherOrderStress.size() == dim * dim * dim, "The referenceHigherOrder stress must have a size of 27" );

        term3 = variableVector( dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        term3[ dim * I + J ] += referenceHigherOrderStress[ dim * dim * I + dim * Q + R ] * invCGamma[ dim * dim * J + dim * Q + R ];
                    }
                }
            }
        }

        return;
    }

    void computeLinearElasticTerm3( const variableVector &invCGamma,
                                        const variableVector &referenceHigherOrderStress, variableVector &term3,
                                        variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress ){
        /*!
         * Compute the value of the third term in the micromorphic linear elasticity formulation.
         * \f$term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR}\f$
         *
         * Also returns the Jacobians
         * \f$\frac{ \partial term3_{IJ} }{ \partial M_{TUV} } = \delta_{IT} C_{JS}^{-1} \Gamma_{SUV}\f$
         * \f$\frac{ \partial term3_{IJ} }{ \partial C_{TW}^{-1} \Gamma_{WUV} = M_{IUV} \delta_{JT}\f$
         *
         * \param &invCGamma: \f$C_{JS}^{-1} \Gamma_{SQR}\f$
         * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * \param &term3: The third term in the linear elastic equation.
         * \param &dTerm3dInvCGamma: The derivative of the third term w.r.t. the product of the inverse RCG tensor and Gamma
         * \param &dTerm3dReferenceHigherOrderStress: The derivative of the third term w.r.t. the higher order stress
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress, term3 ) );

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dTerm3dInvCGamma = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );
        dTerm3dReferenceHigherOrderStress = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int T = 0; T < dim; T++ ){
                    for ( unsigned int U = 0; U < dim; U++ ){
                        for ( unsigned int V = 0; V < dim; V++ ){
                            dTerm3dInvCGamma[ dim * I + J ][ dim * dim * T + dim * U + V ] = referenceHigherOrderStress[ dim * dim * I + dim * U + V ] * eye[ dim * J + T ];
                            dTerm3dReferenceHigherOrderStress[ dim * I + J ][ dim * dim * T + dim * U + V ] = eye[ dim * I + T ] * invCGamma[ dim * dim * J + dim * U + V ];
                        }
                    }
                }
            }
        }
        return;
    }

    void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi ){
        /*!
         * Compute the product \f$C_{IK}^{-1} \Psi_{KJ}\f$
         *
         * \param &invRCG: The inverse of the right cauchy green deformation tensor.
         * \param &Psi: The micro-deformation measure.
         * \param &invRCGPsi: the product.
         */

        //Assume 3d
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( invRCG.size() == dim * dim, "invRCG has an improper dimension" );

        TARDIGRADE_ERROR_TOOLS_CHECK( Psi.size() == dim * dim, "Psi has an improper dimension" );

        invRCGPsi = tardigradeVectorTools::matrixMultiply( invRCG, Psi, dim, dim, dim, dim );

        return;
    }

    void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                               variableMatrix &dInvRCGPsidRCG, variableMatrix &dInvRCGPsidPsi ){
        /*!
         * Compute the product \f$C_{IK}^{-1} \Psi_{KJ}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial C_{KL} } = -C_{IK}^{-1} C_{LO}^{-1} \Psi_{OJ}\f$
         * \f$\frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial \Psi_{KL} } = C_{IK}^{-1} \delta_{JL}\f$
         *
         * \param &invRCG: The inverse of the right cauchy green deformation tensor.
         * \param &Psi: The micro-deformation measure.
         * \param &invRCGPsi: the product.
         * \param &dInvRCGPsidRCG: The Jacobian of the product w.r.t. the right cauchy green
         *     deformation tensor.
         * \param &dInvRCGPsidPsi: The Jacobian of the product w.r.t. the micro-deformation measure.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi ) );

        //Construct the jacobians
        variableVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dInvRCGPsidRCG = variableMatrix( invRCGPsi.size(), variableVector( invRCG.size(), 0 ) );
        dInvRCGPsidPsi = variableMatrix( invRCGPsi.size(), variableVector( Psi.size(), 0 ) );

        for ( unsigned int I = 0; I < 3; I++ ){
            for ( unsigned int J = 0; J < 3; J++ ){
                for ( unsigned int K = 0; K < 3; K++ ){
                    for ( unsigned int L = 0; L < 3; L++ ){
                        dInvRCGPsidRCG[ dim * I + J ][ dim * K + L ] = -invRCG[ dim * I + K ] * invRCGPsi[ dim * L + J ];
                        dInvRCGPsidPsi[ dim * I + J ][ dim * K + L ] = invRCG[ dim * I + K ] * eye[ dim * J + L ];
                    }
                }
            }
        }
        
        return;
    }

    void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma ){
        /*!
         * Compute the product \f$C_{IS}^{-1} \Gamma_{SQR}\f$
         *
         * \param &invRCG: The inverse of the right Cauchy Green deformation tensor.
         * \param &Gamma: The gradient of the micro-deformation deformation tensor.
         * \param &invRCGGamma: The product.
         */

        //Assume 3d
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( invRCG.size() == dim * dim, "invRCG has an improper dimension" );

        TARDIGRADE_ERROR_TOOLS_CHECK( Gamma.size() == dim * dim * dim, "Gamma has an improper dimension" );

        invRCGGamma = variableVector( dim * dim * dim, 0 );
        for ( unsigned int J = 0; J < dim; J++ ){
            for ( unsigned int Q = 0; Q < dim; Q++ ){
                for ( unsigned int R = 0; R < dim; R++ ){
                    for ( unsigned int S = 0; S < dim; S++ ){
                        invRCGGamma[ dim * dim * J + dim * Q + R ] += invRCG[ dim * J + S ] * Gamma[ dim * dim * S + dim * Q + R ];
                    }
                }
            }
        }

        return;
    }

    void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                 variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma ){
        /*!
         * Compute the product \f$C_{IS}^{-1} \Gamma_{SQR}\f$
         *
         * Also compute the Jacobians
         * \f$\frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial C_{TU} } = -C_{JT}^{-1} C_{US}^{-1} \Gamma_{SQR}\f$
         * \f$\frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial \Gamma_{TUV} } = C_{JT}^{-1} \delta_{QU} \delta_{RV}\f$
         *
         * \param &invRCG: The inverse of the right Cauchy Green deformation tensor.
         * \param &Gamma: The gradient of the micro-deformation deformation tensor.
         * \param &invRCGGamma: The product.
         * \param &dInvRCGGammadRCG: The derivative of the inverse RCG gamma product w.r.t. the RCG tensor
         * \param &dInvRCGGammadGamma: The derivative of the inverse RCG gamma product w.r.t. Gamma
         */

        //Assume 3d
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma ) );

        //Assemble jacobians of invCGamma w.r.t. C and Gamma
        variableVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        dInvRCGGammadRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dInvRCGGammadGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int J = 0; J < dim; J++ ){
            for ( unsigned int Q = 0; Q < dim; Q++ ){
                for ( unsigned int R = 0; R < dim; R++ ){
                    for ( unsigned int T = 0; T < dim; T++ ){
                        for ( unsigned int U = 0; U < dim; U++ ){
                            dInvRCGGammadRCG[ dim * dim * J + dim * Q + R ][ dim * T + U ] 
                                = -invRCG[ dim * J + T] * invRCGGamma[ dim * dim * U + dim * Q + R ];
                            for ( unsigned int V = 0; V < dim; V++ ){
                                dInvRCGGammadGamma[ dim * dim * J + dim * Q + R ][ dim * dim * T + dim * U + V]
                                    = invRCG[ dim * J + T ] * eye[ dim * Q + U ] * eye[ dim * R + V ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress ){
        /*!
         * Map the stress measures in the reference configuration to the current configuration.
         *
         * \param &deformationGradient: The deformation gradient between the 
         *     reference configuration and the current configuration.
         * \param &microDeformation: The micro-deformation map between the 
         *     reference configuration and the current configuration.
         * \param &PK2Stress: The Second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher order stress in 
         *     the reference configuration.
         * \param &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
         * \param &microStress: The symmetric micro-stress in the current configuration.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         */

        //Map the PK2 stress to the Cauchy stress
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress ) );

        //Map the symmetric micro stress to the current configuration
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress ) );

        //Map the higher order stress to the current configuration
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                                                 microDeformation, higherOrderStress ) );

        return;
    }

    void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress,
                                         variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                         variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                         variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                         variableMatrix &dHigherOrderStressdReferenceHigherOrderStress ){
        /*!
         * Map the stress measures in the reference configuration to the current configuration.
         *
         * Also computes the Jacobians
         *
         * \param &deformationGradient: The deformation gradient between the 
         *     reference configuration and the current configuration.
         * \param &microDeformation: The micro-deformation map between the 
         *     reference configuration and the current configuration.
         * \param &PK2Stress: The Second Piola-Kirchoff stress.
         * \param &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * \param &referenceHigherOrderStress: The higher order stress in 
         *     the reference configuration.
         * \param &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
         * \param &microStress: The symmetric micro-stress in the current configuration.
         * \param &higherOrderStress: The higher order stress in the current configuration.
         * \param &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the 
         *     deformation gradient.
         * \param &dCauchyStressdPK2Stress: The Jacobian of the Cauchy stress w.r.t. the 
         *     PK2 stress.
         * \param &dMicroStressdF: The Jacobian of the micro stress w.r.t. the 
         *     deformation gradient.
         * \param &dMicroStressdReferenceMicroStress: The Jacobian of the micro-stress 
         *     in the current configuration w.r.t. the micro-stress in the reference configuration.
         * \param &dHigherOrderStressdF: The Jacobian of the higher-order stress w.r.t.
         *     the deformation gradient.
         * \param &dHigherOrderStressdChi: The Jacobian of the higher-order stress 
         *     w.r.t. the micro-deformation.
         * \param &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     higher-order stress w.r.t. the higher order stress in the reference configuration.
         */

        //Map the PK2 stress to the Cauchy stress
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress,
                                                                                         dCauchyStressdPK2Stress, dCauchyStressdF ) );

        //Map the symmetric micro stress to the current configuration
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress,
                                                                                                    dMicroStressdReferenceMicroStress, dMicroStressdF ) );

        //Map the higher order stress to the current configuration
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                                                 microDeformation, higherOrderStress,
                                                                                                 dHigherOrderStressdReferenceHigherOrderStress,
                                                                                                 dHigherOrderStressdF,
                                                                                                 dHigherOrderStressdChi ) );

        return;
    }

    void linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress ){
        /*!
         * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &cauchyStress: The Cauchy stress.
         * \param &microStress: The symmetric micro-stress.
         * \param &higherOrderStress: The higher-order stress.
         */

        variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
        TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                 A, B, C, D,
                                                                 PK2Stress, referenceMicroStress, referenceHigherOrderStress ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                                  referenceMicroStress, referenceHigherOrderStress,
                                                                  cauchyStress, microStress, higherOrderStress ) );

        return;
    }
    
    void linearElasticity(  const variableVector &deformationGradient, const variableVector &microDeformation,
                                const variableVector &gradientMicroDeformation,
                                const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                const parameterVector &D,
                                variableVector &cauchyStress, variableVector &microStress,
                                variableVector &higherOrderStress,
                                variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdChi, variableMatrix &dCauchyStressdGradChi,
                                variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdChi, variableMatrix &dMicroStressdGradChi,
                                variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                variableMatrix &dHigherOrderStressdGradChi ){
        /*!
         * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * Also compute the Jacobians
         *
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro-deformation
         * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * \param &A: The A stiffness matrix.
         * \param &B: The B stiffness matrix.
         * \param &C: The C stiffness matrix.
         * \param &D: The D stiffness matrix.
         * \param &cauchyStress: The Cauchy stress.
         * \param &microStress: The symmetric micro-stress.
         * \param &higherOrderStress: The higher-order stress.
         * \param &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the deformation gradient
         * \param &dCauchyStressdChi: The Jacobian of the Cauchy stress w.r.t. the micro deformation.
         * \param &dCauchyStressdGradChi: The Jacobian of the Cauchy stress w.r.t. the gradient of the 
         *     micro-deformation.
         * \param &dMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient
         * \param &dMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
         * \param &dMicroStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the 
         *     micro-deformation.
         * \param &dHigherOrderStressdF: The Jacobian of the Higher Order stress w.r.t. the deformation gradient
         * \param &dHigherOrderStressdChi: The Jacobian of the Higher Order stress w.r.t. the micro deformation.
         * \param &dHigherOrderStressdGradChi: The Jacobian of the Higher Order stress w.r.t. the gradient of the 
         *     micro-deformation.
         */

        variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;

        variableMatrix dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi;
        variableMatrix dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi;
        variableMatrix dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi;

        TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                 A, B, C, D,
                                                                 PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                 dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi,
                                                                 dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi,
                                                                 dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi ) );

        variableMatrix dCauchyStressdPK2Stress, dMicroStressdReferenceMicroStress, dHigherOrderStressdReferenceHigherOrderStress;

        TARDIGRADE_ERROR_TOOLS_CATCH( mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                                  referenceMicroStress, referenceHigherOrderStress,
                                                                  cauchyStress, microStress, higherOrderStress,
                                                                  dCauchyStressdF, dCauchyStressdPK2Stress,
                                                                  dMicroStressdF, dMicroStressdReferenceMicroStress,
                                                                  dHigherOrderStressdF, dHigherOrderStressdChi,
                                                                  dHigherOrderStressdReferenceHigherOrderStress ) );

        //Assemble the jacobians of the Cauchy stress
        dCauchyStressdF += tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdF );
        dCauchyStressdChi = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdChi );
        dCauchyStressdGradChi = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdGradChi );

        //Assemble the jacobians of the symmetric micro-stress
        dMicroStressdF += tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdF );
        dMicroStressdChi = tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdChi );
        dMicroStressdGradChi = tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdGradChi );

        //Assemble the jacobians of the higher-order stress
        dHigherOrderStressdF += tardigradeVectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                  dReferenceHigherOrderStressdF );
        dHigherOrderStressdGradChi = tardigradeVectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                      dReferenceHigherOrderStressdGradChi );

        return;
    }

    void formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A ){
        /*!
         * Form the isotropic A stiffness tensor.
         * \f$A_{KLMN} = \lambda \delta_{KL} \delta_{MN} + \mu \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )\f$
         *
         * \param &lambda: The micromorphic lambda parameter.
         * \param &mu: The micromorphic mu parameter.
         * \param &A: The isotropic A stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        A = parameterVector( dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        A[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = lambda * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                               + mu * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                      + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return;
    }

    void formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                             const parameterType &nu,  const parameterType &sigma, parameterVector &B ){
        /*!
         * Form the isotropic B stiffness tensor.
         * \f$B_{KLMN} = ( eta - tau ) \delta_{KL} \delta_{MN} + \kappa \delta_{KM} \delta_{LN} + \nu \delta_{KN} \delta_{LM}
         *             - \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right)\f$
         *
         * \param &eta: The micromorphic eta parameter.
         * \param &tau: The micromorphic tau parameter.
         * \param &kappa: The micromorphic kappa parameter.
         * \param &nu: The micromorphic nu parameter.
         * \param &sigma: The micromorphic sigma parameter
         * \param &B: The isotropic B stiffnes tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        B = parameterVector( dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        B[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = ( eta - tau ) * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                               + kappa * eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                               + nu * eye[ dim * K + N ] * eye[ dim * L + M ]
                                                                               - sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                         + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return;
    }

    void formIsotropicC( const parameterVector &taus, parameterVector &C ){
        /*!
         * Form the isotropic C stiffness tensor.
         * \f$\begin{align}
         * C_{KLMNPQ} &= \tau_1 \left( \delta_{KL} \delta_{MN} \delta_{PQ} + \delta_{KQ} \delta_{LM} \delta_{NP} \right)\\
         *            &+ \tau_2 \left( \delta_{KL} \delta_{MP} \delta_{NQ} + \delta_{KM} \delta_{LQ} \delta_{NP} \right)\\
         *            &+ \tau_3 \delta_{KL} \delta_{MQ} \delta_{NP}\\
         *            &+ \tau_4 \delta_{KN} \delta_{LM} \delta_{PQ}\\
         *            &+ \tau_5 \left( \delta_{KM} \delta_{LN} \delta_{PQ} + \delta_{KP} \delta_{LM} \delta_{NQ} )\\
         *            &+ \tau_6 \delta_{KM} \delta_{LP} \delta_{NQ}\\
         *            &+ \tau_7 \delta_{KN} \delta_{LP} \delta_{MQ}\\
         *            &+ \tau_8 \left( \delta_{KP} \delta_{LQ} \delta_{MN} + \delta_{KQ} \delta_{LN} \delta_{MP} )\\
         *            &+ \tau_9 \delta_{KN} \delta_{LQ} \delta_{MP}\\
         *            &+ \tau_{10} \delta_{KP} \delta_{LN} \delta_{MQ}\\
         *            &+ \tau_{11} \delta_{KQ} \delta_{LP} \delta_{MN}
         * \end{align}\f$
         *
         * \param &taus: The moduli (11 independent terms)
         * \param &C: The isotropic C stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        TARDIGRADE_ERROR_TOOLS_CHECK( taus.size() == 11, "11 moduli required to form C" );

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        C = parameterVector( dim * dim * dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int P = 0; P < dim; P++ ){
                            for ( unsigned int Q = 0; Q < dim; Q++ ){
                                C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M 
                                 + dim * dim * N + dim * P + Q ] = taus[0] * ( eye[ dim * K + L ] * eye[ dim * M + N ] * eye[ dim * P + Q ]
                                                                             + eye[ dim * K + Q ] * eye[ dim * L + M ] * eye[ dim * N + P ] )
                                                                 + taus[1] * ( eye[ dim * K + L ] * eye[ dim * M + P ] * eye[ dim * N + Q ]
                                                                             + eye[ dim * K + M ] * eye[ dim * L + Q ] * eye[ dim * N + P ] )
                                                                 + taus[2] * eye[ dim * K + L ] * eye[ dim * M + Q ] * eye[ dim * N + P]
                                                                 + taus[3] * eye[ dim * K + N ] * eye[ dim * L + M ] * eye[ dim * P + Q]
                                                                 + taus[4] * ( eye[ dim * K + M ] * eye[ dim * L + N ] * eye[ dim * P + Q ]
                                                                             + eye[ dim * K + P ] * eye[ dim * L + M ] * eye[ dim * N + Q ] )
                                                                 + taus[5] * eye[ dim * K + M ] * eye[ dim * L + P ] * eye[ dim * N + Q ]
                                                                 + taus[6] * eye[ dim * K + N ] * eye[ dim * L + P ] * eye[ dim * M + Q ]
                                                                 + taus[7] * ( eye[ dim * K + P ] * eye[ dim * L + Q ] * eye[ dim * M + N ]
                                                                             + eye[ dim * K + Q ] * eye[ dim * L + N ] * eye[ dim * M + P ] )
                                                                 + taus[8] * eye[ dim * K + N ] * eye[ dim * L + Q ] * eye[ dim * M + P ]
                                                                 + taus[9] * eye[ dim * K + P ] * eye[ dim * L + N ] * eye[ dim * M + Q ]
                                                                 + taus[10] * eye[ dim * K + Q ] * eye[ dim * L + P ] * eye[ dim * M + N ];
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    void formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D ) {
        /*!
         * Form the isotropic tensor D.
         * \f$D_{KLMN} = \tau \delta_{KL} \delta_{MN} + \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )\f$
         *
         * \param &tau: The micromorphic tau parameter.
         * \param &sigma: The micromorphic sigma parameter.
         * \param &D: The D stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        tardigradeVectorTools::eye( eye );

        D = parameterVector( dim * dim * dim * dim, 0 );
        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        D[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = tau * eye[ dim * K + L ] * eye[ dim * M + N ]
                            + sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ] + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return;
    }

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation ){
        /*!
         * Assemble the fundamental deformation meaures from the degrees of freedom.
         *
         * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
         * \param &phi: The micro displacement.
         * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro deformation
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         */


        //Extract the degrees of freedom
        variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                     grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                     grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                     grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                     grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                     grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                     grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                     grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                     grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation ) );

        return;
    }

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation, variableMatrix &dFdGradU,
                                                     variableMatrix &dChidPhi, variableMatrix &dGradChidGradPhi ){
        /*!
         * Assemble the fundamental deformation meaures from the degrees of freedom.
         *
         * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
         * \param &phi: The micro displacement.
         * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro deformation
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         * \param &dFdGradU: The Jacobian of the deformation gradient w.r.t. the gradient of the displacement
         * \param &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
         * \param &dGradChidGradPhi: The Jacobian of the gradient of the micro deformation w.r.t.
         *      the gradient of the micro displacement
         */

        const unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;
        const unsigned int tot_dim = sot_dim * dim;

        //Extract the degrees of freedom
        variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                     grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                     grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                     grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                     grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                     grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                     grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                     grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                     grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };

        variableVector _dFdGradU, _dChidPhi, _dGradChidGradPhi;

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, _dFdGradU ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, _dChidPhi ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
                                                                                                     _dGradChidGradPhi ) );

        dFdGradU         = tardigradeVectorTools::inflate( _dFdGradU,         sot_dim, sot_dim );
        dChidPhi         = tardigradeVectorTools::inflate( _dChidPhi,         sot_dim, sot_dim );
        dGradChidGradPhi = tardigradeVectorTools::inflate( _dGradChidGradPhi, tot_dim, tot_dim );

        return;
    }

    void extractMaterialParameters( const std::vector< double > &fparams,
                                        parameterVector &Amatrix, parameterVector &Bmatrix,
                                        parameterVector &Cmatrix, parameterVector &Dmatrix ){
        /*!
         * Extract the parameters from the parameter vector
         *
         * \param &fparams: The incoming parameter vector
         * \param &Amatrix: The A stiffness matrix.
         * \param &Bmatrix: The B stiffness matrix.
         * \param &Cmatrix: The C stiffness matrix.
         * \param &Dmatrix: The D stiffness matrix.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( fparams.size() != 0, "The material parameters vector has a length of 0" );

        unsigned int start = 0;
        unsigned int span;

        std::vector< parameterVector > outputs( 4 );

        //Extract the material parameters
        for ( unsigned int i = 0; i < outputs.size(); i++ ){
            span = ( unsigned int )std::floor( fparams[ start ]  + 0.5 ); //Extract the span of the parameter set

            TARDIGRADE_ERROR_TOOLS_CHECK( fparams.size() >= start + 1 + span, "fparams is not long enough to contain all of the required parameters:\n    filling variable " + std::to_string( i )
                                                                            + "\n    size =          "  + std::to_string( fparams.size() )
                                                                            + "\n    required size = "  + std::to_string( start + 1 + span ) );

            outputs[ i ] = parameterVector( fparams.begin() + start + 1, fparams.begin() + start + 1 + span );

            start = start + 1 + span;
        }

        //Form the stiffness tensors
        TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 0 ].size( ) == 2, "Unrecognized number of parameters ( " + std::to_string( outputs[ 0 ].size() ) + " ) for the A stiffness tensor" );
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicLinearElasticity::formIsotropicA( outputs[ 0 ][ 0 ], outputs[ 0 ][ 1 ], Amatrix ) );

        TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 1 ].size() == 5, "Unrecognized number of parameters ( " + std::to_string( outputs[ 1 ].size() ) + " ) for the B stiffness tensor" );
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicLinearElasticity::formIsotropicB( outputs[ 1 ][ 0 ], outputs[ 1 ][ 1 ], outputs[ 1 ][ 2 ],
                                                                                              outputs[ 1 ][ 3 ], outputs[ 1 ][ 4 ], Bmatrix ) );

        TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 2 ].size() == 11, "Unrecognized number of parameters ( " + std::to_string( outputs[ 1 ].size() ) + " ) for the B stiffness tensor" );
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicLinearElasticity::formIsotropicC( outputs[ 2 ], Cmatrix ) );

        TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 3 ].size() == 2, "Unrecognized number of parameters ( " + std::to_string( outputs[ 3 ].size() ) + " ) for the D stiffness tensor" );
        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicLinearElasticity::formIsotropicD( outputs[ 3 ][ 0 ], outputs[ 3 ][ 1 ], Dmatrix ) );

        return;
    }


    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ), 
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::string &output_message
                      ){
        /*!
         * Evaluate the linear-elastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacement values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables ( unused )
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback. ( unused )
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Re-direct the output to a buffer
        std::stringbuf buffer;
        cerr_redirect rd( &buffer );

        /*=============================
        | Extract the incoming values |
        ==============================*/

        //Extract the parameters
        parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;

        try{
            extractMaterialParameters( fparams, Amatrix, Bmatrix, Cmatrix, Dmatrix );
        }
        catch( std::exception &e ){
            tardigradeErrorTools::captureNestedExceptions( e, output_message );
            return 2;
        }

        /*===============================================
        | Assemble the fundamental deformation measures |
        ================================================*/

        //Compute the fundamental deformation measures from the degrees of freedom

        variableVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

        try{
            assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                    currentDeformationGradient, currentMicroDeformation,
                                                    currentGradientMicroDeformation );
        }
        catch( std::exception &e ){

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            return 2;
        }

        /*===============================
        | Compute the new stress values |
        ===============================*/

        //Compute the new stress values
        variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

        try{

            tardigradeMicromorphicLinearElasticity::linearElasticityReference( currentDeformationGradient,
                                                                               currentMicroDeformation,
                                                                               currentGradientMicroDeformation,
                                                                               Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                               PK2, SIGMA, M );

        }
        catch( std::exception &e ){

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            return 2;

        }

        //No errors in calculation.
        return 0;
    }

    int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ), 
                        const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                        const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                        const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                        std::vector< double > &SDVS,
                        const std::vector< double > &current_ADD_DOF,
                        const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                        const std::vector< double > &previous_ADD_DOF,
                        const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                        std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                        std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                        std::vector< std::vector< double > > &DPK2Dgrad_phi, 
                        std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                        std::vector< std::vector< double > > &DSIGMADgrad_phi,
                        std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                        std::vector< std::vector< double > > &DMDgrad_phi,
                        std::vector< std::vector< double > > &ADD_TERMS,
                        std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                        std::string &output_message
                      ){
        /*!
         * Evaluate the linear-elastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacement values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         * \param &PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the 
         *     gradient of macro displacement.
         * \param &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * \param &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * \param &DSIGMADgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * \param &DSIGMADphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * \param &DSIGMADgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * \param &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * \param &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * \param &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation ( unused )
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        //Re-direct the output to a buffer
        std::stringbuf buffer;
        cerr_redirect rd( &buffer );

        /*=============================
        | Extract the incoming values |
        ==============================*/

        //Extract the parameters
        parameterVector Amatrix, Bmatrix, Cmatrix, Dmatrix;

        try{
            extractMaterialParameters( fparams, Amatrix, Bmatrix, Cmatrix, Dmatrix );
        }
        catch( std::exception &e ){
            tardigradeErrorTools::captureNestedExceptions( e, output_message );
            return 2;
        }

        /*===============================================
        | Assemble the fundamental deformation measures |
        ================================================*/

        //Compute the fundamental deformation measures from the degrees of freedom

        variableVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;
        variableMatrix dDeformationGradientdGradU, dMicroDeformationdPhi, dGradientMicroDeformationdGradPhi;

        try{
            assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                    currentDeformationGradient, currentMicroDeformation,
                                                    currentGradientMicroDeformation, dDeformationGradientdGradU,
                                                    dMicroDeformationdPhi, dGradientMicroDeformationdGradPhi );
        }
        catch( std::exception &e ){
            tardigradeErrorTools::captureNestedExceptions( e, output_message );
            return 2;
        }

        /*===============================
        | Compute the new stress values |
        ===============================*/

        //Compute the new stress values
        variableVector currentPK2Stress, currentReferenceMicroStress, currentReferenceHigherOrderStress;

        variableMatrix dPK2dDeformationGradient, dPK2dMicroDeformation, dPK2dGradientMicroDeformation,
                       dSIGMAdDeformationGradient, dSIGMAdMicroDeformation, dSIGMAdGradientMicroDeformation,
                       dMdDeformationGradient, dMdGradientMicroDeformation;

        try{
            tardigradeMicromorphicLinearElasticity::linearElasticityReference( currentDeformationGradient,
                                                                     currentMicroDeformation,
                                                                     currentGradientMicroDeformation,
                                                                     Amatrix, Bmatrix, Cmatrix, Dmatrix,
                                                                     PK2, SIGMA, M,
                                                                     dPK2dDeformationGradient, dPK2dMicroDeformation,
                                                                     dPK2dGradientMicroDeformation,
                                                                     dSIGMAdDeformationGradient, dSIGMAdMicroDeformation,
                                                                     dSIGMAdGradientMicroDeformation,
                                                                     dMdDeformationGradient, dMdGradientMicroDeformation );
        }
        catch( std::exception &e ){
            tardigradeErrorTools::captureNestedExceptions( e, output_message );
            return 2;
        }

        /*=======================
        | Assemble the Jacobian |
        =======================*/

        DPK2Dgrad_u     = tardigradeVectorTools::dot( dPK2dDeformationGradient, dDeformationGradientdGradU );
        DPK2Dphi        = tardigradeVectorTools::dot( dPK2dMicroDeformation, dMicroDeformationdPhi );
        DPK2Dgrad_phi   = tardigradeVectorTools::dot( dPK2dGradientMicroDeformation, dGradientMicroDeformationdGradPhi );

        DSIGMADgrad_u   = tardigradeVectorTools::dot( dSIGMAdDeformationGradient, dDeformationGradientdGradU );
        DSIGMADphi      = tardigradeVectorTools::dot( dSIGMAdMicroDeformation, dMicroDeformationdPhi );
        DSIGMADgrad_phi = tardigradeVectorTools::dot( dSIGMAdGradientMicroDeformation, dGradientMicroDeformationdGradPhi );

        DMDgrad_u       = tardigradeVectorTools::dot( dMdDeformationGradient, dDeformationGradientdGradU );
        DMDphi          = variableMatrix( 27, variableVector( 9, 0 ) );
        DMDgrad_phi     = tardigradeVectorTools::dot( dMdGradientMicroDeformation, dGradientMicroDeformationdGradPhi );

        //No errors in calculation.
        return 0;
    }

    //! Define the hydra version of the micromorphic linear elasticity model
    class hydraMicromorphicLinearElasticity : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            //! The elasticity residual class
            tardigradeHydra::micromorphicLinearElasticity::residual elasticity;

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{
                /*!
                 * Set the vector of residual classes (in this case, only elasticity)
                 */

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                elasticity = tardigradeHydra::micromorphicLinearElasticity::residual( this, *getConfigurationUnknownCount( ), *getParameters( ) );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    void generate_input_variable_string( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                         const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                         const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                         const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                         std::vector< double > &SDVS,
                                         const std::vector< double > &current_ADD_DOF,
                                         const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                         const std::vector< double > &previous_ADD_DOF,
                                         const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                         std::string &input_variables ){
        /*
         * Summarize the input variables in string form for debugging
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacment values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &input_variables: The input variables in string form
         */

        input_variables = "";

        input_variables += "time:\n";
        for ( auto t = time.begin( ); t != time.end( ); t++ ){ input_variables += " " + std::to_string( *t ) + ","; }
        input_variables += "\nfparams:\n";
        for ( auto f = fparams.begin( ); f != fparams.end( ); f++ ){ input_variables += " " + std::to_string( *f ) + ","; }
        input_variables += "\ncurrent_grad_u:\n";
        for ( unsigned int i = 0; i < 3; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ input_variables += " " + std::to_string( current_grad_u[ i ][ j ] ) + ","; } input_variables += "\n"; }
        input_variables += "\ncurrent_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ input_variables += " " + std::to_string( current_phi[ i ] ) + ","; }
        input_variables += "\ncurrent_grad_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ input_variables += " " + std::to_string( current_grad_phi[ i ][ j ] ) + ","; } input_variables += "\n"; }
        input_variables += "\nprevious_grad_u:\n";
        for ( unsigned int i = 0; i < 3; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ input_variables += " " + std::to_string( previous_grad_u[ i ][ j ] ) + ","; } input_variables += "\n"; }
        input_variables += "\nprevious_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ input_variables += " " + std::to_string( previous_phi[ i ] ) + ","; }
        input_variables += "\nprevious_grad_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ input_variables += " " + std::to_string( previous_grad_phi[ i ][ j ] ) + ","; } input_variables += "\n"; }
        input_variables += "\nSDVS:\n";
        for ( auto s = SDVS.begin( ); s != SDVS.end( ); s++ ){ input_variables += " " + std::to_string( *s ) + ","; }
        input_variables += "\ncurrent_ADD_DOF:\n";
        for ( auto a = current_ADD_DOF.begin( ); a != current_ADD_DOF.end( ); a++ ){ input_variables += " " + std::to_string( *a ) + ","; }
        input_variables += "\ncurrent_ADD_grad_DOF:\n";
        for ( auto a = current_ADD_grad_DOF.begin( ); a != current_ADD_grad_DOF.end( ); a++ ){ for ( auto g = a->begin( ); g != a->end( ); g++ ){ input_variables += " " + std::to_string( *g ) + ","; } input_variables += "\n"; }
        input_variables += "\nprevious_ADD_DOF:\n";
        for ( auto a = previous_ADD_DOF.begin( ); a != previous_ADD_DOF.end( ); a++ ){ input_variables += " " + std::to_string( *a ) + ","; }
        input_variables += "\nprevious_ADD_grad_DOF:\n";
        for ( auto a = previous_ADD_grad_DOF.begin( ); a != previous_ADD_grad_DOF.end( ); a++ ){ for ( auto g = a->begin( ); g != a->end( ); g++ ){ input_variables += " " + std::to_string( *g ) + ","; } input_variables += "\n"; }

        return;

    }

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ), 
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::string &output_message
                            ){
        /*!
         * Evaluate the linear-elastic constitutive model using the hydra computational framework.
         * Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacement values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables ( unused )
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         * \param &PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback. ( unused )
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        variableType temperature         = 293.15; // Tardigrade doesn't have temperature for micromorphic currently so we're hardcoding these
        variableType previousTemperature = 293.15;

        variableVector currentDeformationGradient,  currentMicroDeformation,  currentGradientMicroDeformation;
        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

        try{

            /*===============================================
            | Assemble the fundamental deformation measures |
            ================================================*/

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation )
            )

            hydraMicromorphicLinearElasticity hydra( time[ 0 ], time[ 1 ],
                                                     temperature,                     previousTemperature,
                                                     currentDeformationGradient,      previousDeformationGradient,
                                                     currentMicroDeformation,         previousMicroDeformation,
                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                     { }, { },
                                                     SDVS, fparams, 1, 0, 3, 45, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

            // Compute the stress
            hydra.evaluate( );

            PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                    hydra.getUnknownVector( )->begin( ) +  9 );

            SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                    hydra.getUnknownVector( )->begin( ) + 18 );

            M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                    hydra.getUnknownVector( )->begin( ) + 45 );

        }
        catch( tardigradeHydra::convergence_error &e ){

            //Convergence error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            output_message = "INPUT PARAMETERS FOLLOW:\n" + input_variables;

            return 1;

        }
        catch( std::exception &e ){

            //Fatal error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables;

#ifdef TARDIGRADE_FATAL_AS_CONVERGENCE
            return 1;
#else
            return 2;
#endif

        }

        //No errors in calculation.
        return 0;
    }

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ), 
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                              std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                              std::vector< std::vector< double > > &DPK2Dgrad_phi, 
                              std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                              std::vector< std::vector< double > > &DSIGMADgrad_phi,
                              std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                              std::vector< std::vector< double > > &DMDgrad_phi,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                              std::string &output_message
                            ){
        /*!
         * Evaluate the linear-elastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacement values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the 
         *     gradient of macro displacement.
         * \param &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * \param &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * \param &DSIGMADgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * \param &DSIGMADphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * \param &DSIGMADgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * \param &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * \param &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * \param &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation ( unused )
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        variableType temperature         = 293.15; // Tardigrade doesn't have temperature for micromorphic currently so we're hardcoding these
        variableType previousTemperature = 293.15;

        variableVector currentDeformationGradient,  currentMicroDeformation,  currentGradientMicroDeformation;
        variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;
        variableMatrix previousdFdGradU, previousdChidPhi, previousdGradChidGradPhi;

        try{

            /*===============================================
            | Assemble the fundamental deformation measures |
            ================================================*/

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation,
                                                        dFdGradU, dChidPhi, dGradChidGradPhi )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation,
                                                        previousdFdGradU, previousdChidPhi, previousdGradChidGradPhi )
            )

            hydraMicromorphicLinearElasticity hydra( time[ 0 ], time[ 1 ],
                                                     temperature,                     previousTemperature,
                                                     currentDeformationGradient,      previousDeformationGradient,
                                                     currentMicroDeformation,         previousMicroDeformation,
                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                     { }, { },
                                                     SDVS, fparams, 1, 0, 3, 45, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

            // Compute the stress
            hydra.evaluate( );

            PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                    hydra.getUnknownVector( )->begin( ) +  9 );

            SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                    hydra.getUnknownVector( )->begin( ) + 18 );

            M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                    hydra.getUnknownVector( )->begin( ) + 45 );

            // Compute the consistent tangents
            hydra.computeTangents( );
            const variableVector *dXdD = hydra.getFlatdXdD( );

            unsigned int numConfigurationUnknowns = *hydra.getConfigurationUnknownCount( );

            DPK2Dgrad_u     = variableMatrix(  9, variableVector( 9, 0 ) );

            DSIGMADgrad_u   = variableMatrix(  9, variableVector( 9, 0 ) );

            DMDgrad_u       = variableMatrix( 27, variableVector( 9, 0 ) );

            DPK2Dphi        = variableMatrix(  9, variableVector( 9, 0 ) );

            DSIGMADphi      = variableMatrix(  9, variableVector( 9, 0 ) );

            DMDphi          = variableMatrix( 27, variableVector( 9, 0 ) );

            DPK2Dgrad_phi   = variableMatrix(  9, variableVector( 27, 0 ) );

            DSIGMADgrad_phi = variableMatrix(  9, variableVector( 27, 0 ) );

            DMDgrad_phi     = variableMatrix( 27, variableVector( 27, 0 ) );

            for ( unsigned int i = 0; i < 9; i++ ){

                for ( unsigned int j = 0; j < 9; j++ ){

                    for ( unsigned int k = 0; k < 9; k++ ){

                        DPK2Dgrad_u[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 0 + k ] * dFdGradU[ k ][ j ];

                        DPK2Dphi[ i ][ j ]      += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 9 + k ] * dChidPhi[ k ][ j ];

                        DSIGMADgrad_u[ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 0 + k ] * dFdGradU[ k ][ j ];

                        DSIGMADphi[ i ][ j ]    += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 9 + k ] * dChidPhi[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    for ( unsigned int k = 0; k < 27; k++ ){

                        DPK2Dgrad_phi[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                        DSIGMADgrad_phi[ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    }

                }

            }

            for ( unsigned int i = 0; i < 27; i++ ){

                for ( unsigned int j = 0; j < 9; j++ ){

                    for ( unsigned int k = 0; k < 9; k++ ){

                        DMDgrad_u[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 0 + k ] * dFdGradU[ k ][ j ];

                        DMDphi[ i ][ j ]      += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 9 + k ] * dChidPhi[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < 27; j++ ){

                    for ( unsigned int k = 0; k < 27; k++ ){

                        DMDgrad_phi[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    }

                }

            }

        }
        catch( tardigradeHydra::convergence_error &e ){

            //Convergence error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            output_message = "INPUT PARAMETERS FOLLOW:\n" + input_variables;

            return 1;

        }
        catch( std::exception &e ){

            //Fatal error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables;

#ifdef TARDIGRADE_FATAL_AS_CONVERGENCE
            return 1;
#else
            return 2;
#endif
        }

        //No errors in calculation.
        return 0;

    }

}

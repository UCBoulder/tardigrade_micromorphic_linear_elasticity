/*
 * tardigrade_micromorphic_linear_elasticity.h
 *
 * An implimentation of linear elasticity in the micromorphic context.
 *
 * Based around a quadratic form of the Helmholtz free energy:
 * \rho \psi = \frac{1}{2} E_{IJ} A_{IJKL} E_{KL} + \frac{1}{2} \mathcal{E}_{IJ} B_{IJKL} \mathcal{E}_{KL} 
 *           + \frac{1}{2} \Gamma_{IJK} C_{IJKLMN} \Gamma_{LMN} + E_{IJ} D_{IJKL} \mathcal{E}_{KL}
 */

#ifndef TARDIGRADE_MICROMORPHIC_LINEAR_ELASTICITY_H
#define TARDIGRADE_MICROMORPHIC_LINEAR_ELASTICITY_H

#include<tardigrade_error_tools.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
//#include<micromorphic_material_library.h>

namespace tardigradeMicromorphicLinearElasticity{

    typedef tardigradeMicromorphicTools::variableType variableType;
    typedef tardigradeMicromorphicTools::variableVector variableVector;
    typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

    typedef tardigradeMicromorphicTools::parameterType parameterType;
    typedef tardigradeMicromorphicTools::parameterVector parameterVector;
    typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

    typedef tardigradeMicromorphicTools::constantType constantType;
    typedef tardigradeMicromorphicTools::constantVector constantVector;
    typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

    //! Re-direct std::cout to catch error messages
    struct cout_redirect{
        //! Re-direct the output of std::cout to the buffer
        cout_redirect( std::streambuf * new_buffer)
            : old( std::cout.rdbuf( new_buffer ) )
        { }

        ~cout_redirect( ) {
            std::cout.rdbuf( old );
        }

        private:
            //! A pointer to the buffer for std::cout
            std::streambuf * old;
    };

    //! Re-direct std::cerr to catch error messages
    struct cerr_redirect{
        //! Re-direct the output of std::cerr to the buffer
        cerr_redirect( std::streambuf * new_buffer)
            : old( std::cerr.rdbuf( new_buffer ) )
        { }

        ~cerr_redirect( ) {
            std::cerr.rdbuf( old );
        }

        private:
            //! A pointer to the buffer for std::cerr
            std::streambuf * old;
    };

    void linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress );

    void linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress,
                               variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdChi, variableMatrix &dCauchyStressdGradChi,
                               variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdChi, variableMatrix &dMicroStressdGradChi,
                               variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                               variableMatrix &dHigherOrderStressdGradChi );

    void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress );

    void linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdChi, variableMatrix &dPK2StressdGradChi,
                                        variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdChi,
                                        variableMatrix &dReferenceMicroStressdGradChi, variableMatrix &dMdF, variableMatrix &dMdGradChi );

    void linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                       const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress );

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
                                                       variableMatrix &dMdGamma );

    void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress );

    void mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress,
                                         variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                         variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                         variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                         variableMatrix &dHigherOrderStressdReferenceHigherOrderStress );

    void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma );

    void computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                         variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidChi, 
                                         variableMatrix &dGammadF, variableMatrix &dGammadGradChi );

    void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const parameterVector &A, const parameterVector &D, variableVector &term1 );

    void computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const parameterVector &A, const parameterVector &D, variableVector &term1,
                                        variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain );

    void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &incCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2 );

    void computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                        variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi );

    void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                variableVector &referenceHigherOrderStress );

    void computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress, 
                                                variableMatrix &dHigherOrderStressdGamma );

    void computeLinearElasticTerm3( const variableVector &invCGamma,    
                                        const variableVector &referenceHigherOrderStress, variableVector &term3 );

    void computeLinearElasticTerm3( const variableVector &invCGamma,
                                        const variableVector &referenceHigherOrderStress, variableVector &term3,
                                        variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress );

    void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi );

    void computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                               variableMatrix &dInvRCGPsidRGG, variableMatrix &dInvRCGPsidPsi );

    void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma );

    void computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                 variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma );

    void formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A );

    void formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                             const parameterType &nu,  const parameterType &sigma, parameterVector &B );

    void formIsotropicC( const parameterVector &taus, parameterVector &C );

    void formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D );

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation );

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation,
                                                     variableMatrix &dFdGradU, variableMatrix &dChidPhi,
                                                     variableMatrix &dGradChidGradPhi );

    void extractMaterialParameters( const std::vector< double > &fparams,
                                        parameterVector &Amatrix, parameterVector &Bmatrix,
                                        parameterVector &Cmatrix, parameterVector &Dmatrix );

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
                      );

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
                      );

    void generate_input_variable_string( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                         const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                         const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                         const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                         std::vector< double > &SDVS,
                                         const std::vector< double > &current_ADD_DOF,
                                         const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                         const std::vector< double > &previous_ADD_DOF,
                                         const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                         std::string &input_variables );

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
                            );

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
                            );

}

#endif

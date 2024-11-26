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

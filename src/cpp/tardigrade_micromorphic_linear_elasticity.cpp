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

// Tests for tardigrade_constitutive_tools

#include <tardigrade_micromorphic_linear_elasticity.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_linear_elasticity
#include <boost/test/included/unit_test.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

typedef tardigradeMicromorphicTools::constantType   constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType   parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType   variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

struct cout_redirect {
    cout_redirect(std::streambuf *new_buffer) : old(std::cout.rdbuf(new_buffer)) {}

    ~cout_redirect() { std::cout.rdbuf(old); }

   private:
    std::streambuf *old;
};

struct cerr_redirect {
    cerr_redirect(std::streambuf *new_buffer) : old(std::cerr.rdbuf(new_buffer)) {}

    ~cerr_redirect() { std::cerr.rdbuf(old); }

   private:
    std::streambuf *old;
};

BOOST_AUTO_TEST_CASE(testAssembleFundamentalDeformationMeasures, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     */

    double grad_u[3][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    double phi[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    double grad_phi[9][3] = {
        {1,  2,  3 },
        {4,  5,  6 },
        {7,  8,  9 },
        {10, 11, 12},
        {13, 14, 15},
        {16, 17, 18},
        {19, 20, 21},
        {22, 23, 24},
        {25, 26, 27}
    };

    variableVector answerDeformationGradient = {2, 2, 3, 4, 6, 6, 7, 8, 10};

    variableVector answerMicroDeformation = {2, 2, 3, 4, 6, 6, 7, 8, 10};

    variableVector answerGradientMicroDeformation = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
                                                     15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};

    variableVector resultF, resultChi, resultGradChi;

    tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, phi, grad_phi, resultF,
                                                                                   resultChi, resultGradChi);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultF, answerDeformationGradient));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultChi, answerMicroDeformation));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultGradChi, answerGradientMicroDeformation));

    // Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, phi, grad_phi, resultFJ,
                                                                                   resultChiJ, resultGradChiJ, dFdGradU,
                                                                                   dChidPhi, dGradChidGradPhi);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultFJ, answerDeformationGradient));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultChiJ, answerMicroDeformation));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(resultGradChiJ, answerGradientMicroDeformation));

    // Test the jacobians w.r.t. the gradient of the displacement
    constantType eps = 1e-6;
    for (unsigned int i = 0; i < 9; i++) {
        constantMatrix delta(3, constantVector(3, 0));
        unsigned int   ii, ij;
        ii            = (int)(i / 3);
        ij            = i % 3;
        delta[ii][ij] = eps * fabs(grad_u[ii][ij]) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[3][3] = {
            {grad_u[0][0] + delta[0][0], grad_u[0][1] + delta[0][1], grad_u[0][2] + delta[0][2]},
            {grad_u[1][0] + delta[1][0], grad_u[1][1] + delta[1][1], grad_u[1][2] + delta[1][2]},
            {grad_u[2][0] + delta[2][0], grad_u[2][1] + delta[2][1], grad_u[2][2] + delta[2][2]}
        };

        double negative_perturb[3][3] = {
            {grad_u[0][0] - delta[0][0], grad_u[0][1] - delta[0][1], grad_u[0][2] - delta[0][2]},
            {grad_u[1][0] - delta[1][0], grad_u[1][1] - delta[1][1], grad_u[1][2] - delta[1][2]},
            {grad_u[2][0] - delta[2][0], grad_u[2][1] - delta[2][1], grad_u[2][2] - delta[2][2]}
        };

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(positive_perturb, phi, grad_phi,
                                                                                       FP, chiP, gradChiP);

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(negative_perturb, phi, grad_phi,
                                                                                       FM, chiM, gradChiM);

        variableVector gradCol = (FP - FM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dFdGradU[j][i]);
        }

        gradCol = (chiP - chiM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == 0.);
        }

        gradCol = (gradChiP - gradChiM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == 0.);
        }
    }

    for (unsigned int i = 0; i < 9; i++) {
        constantVector delta(9, 0);

        delta[i] = eps * fabs(phi[i]) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[9] = {phi[0] + delta[0], phi[1] + delta[1], phi[2] + delta[2],
                                      phi[3] + delta[3], phi[4] + delta[4], phi[5] + delta[5],
                                      phi[6] + delta[6], phi[7] + delta[7], phi[8] + delta[8]};

        double negative_perturb[9] = {phi[0] - delta[0], phi[1] - delta[1], phi[2] - delta[2],
                                      phi[3] - delta[3], phi[4] - delta[4], phi[5] - delta[5],
                                      phi[6] - delta[6], phi[7] - delta[7], phi[8] - delta[8]};

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, positive_perturb,
                                                                                       grad_phi, FP, chiP, gradChiP);

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, negative_perturb,
                                                                                       grad_phi, FM, chiM, gradChiM);

        variableVector gradCol = (FP - FM) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(gradCol[j], 0.));
        }

        gradCol = (chiP - chiM) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(gradCol[j], dChidPhi[j][i]));
        }

        gradCol = (gradChiP - gradChiM) / (2 * delta[i]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(gradCol[j], 0.));
        }
    }

    for (unsigned int i = 0; i < 27; i++) {
        constantMatrix delta(9, constantVector(3, 0));
        unsigned int   ii, ij;
        ii            = (int)(i / 3);
        ij            = i % 3;
        delta[ii][ij] = eps * fabs(grad_phi[ii][ij]) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[9][3] = {
            {grad_phi[0][0] + delta[0][0], grad_phi[0][1] + delta[0][1], grad_phi[0][2] + delta[0][2]},
            {grad_phi[1][0] + delta[1][0], grad_phi[1][1] + delta[1][1], grad_phi[1][2] + delta[1][2]},
            {grad_phi[2][0] + delta[2][0], grad_phi[2][1] + delta[2][1], grad_phi[2][2] + delta[2][2]},
            {grad_phi[3][0] + delta[3][0], grad_phi[3][1] + delta[3][1], grad_phi[3][2] + delta[3][2]},
            {grad_phi[4][0] + delta[4][0], grad_phi[4][1] + delta[4][1], grad_phi[4][2] + delta[4][2]},
            {grad_phi[5][0] + delta[5][0], grad_phi[5][1] + delta[5][1], grad_phi[5][2] + delta[5][2]},
            {grad_phi[6][0] + delta[6][0], grad_phi[6][1] + delta[6][1], grad_phi[6][2] + delta[6][2]},
            {grad_phi[7][0] + delta[7][0], grad_phi[7][1] + delta[7][1], grad_phi[7][2] + delta[7][2]},
            {grad_phi[8][0] + delta[8][0], grad_phi[8][1] + delta[8][1], grad_phi[8][2] + delta[8][2]}
        };

        double negative_perturb[9][3] = {
            {grad_phi[0][0] - delta[0][0], grad_phi[0][1] - delta[0][1], grad_phi[0][2] - delta[0][2]},
            {grad_phi[1][0] - delta[1][0], grad_phi[1][1] - delta[1][1], grad_phi[1][2] - delta[1][2]},
            {grad_phi[2][0] - delta[2][0], grad_phi[2][1] - delta[2][1], grad_phi[2][2] - delta[2][2]},
            {grad_phi[3][0] - delta[3][0], grad_phi[3][1] - delta[3][1], grad_phi[3][2] - delta[3][2]},
            {grad_phi[4][0] - delta[4][0], grad_phi[4][1] - delta[4][1], grad_phi[4][2] - delta[4][2]},
            {grad_phi[5][0] - delta[5][0], grad_phi[5][1] - delta[5][1], grad_phi[5][2] - delta[5][2]},
            {grad_phi[6][0] - delta[6][0], grad_phi[6][1] - delta[6][1], grad_phi[6][2] - delta[6][2]},
            {grad_phi[7][0] - delta[7][0], grad_phi[7][1] - delta[7][1], grad_phi[7][2] - delta[7][2]},
            {grad_phi[8][0] - delta[8][0], grad_phi[8][1] - delta[8][1], grad_phi[8][2] - delta[8][2]}
        };

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, phi, positive_perturb,
                                                                                       FP, chiP, gradChiP);

        tardigradeMicromorphicLinearElasticity::assembleFundamentalDeformationMeasures(grad_u, phi, negative_perturb,
                                                                                       FM, chiM, gradChiM);

        variableVector gradCol = (FP - FM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == 0.);
        }

        gradCol = (chiP - chiM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == 0.);
        }

        gradCol = (gradChiP - gradChiM) / (2 * delta[ii][ij]);

        for (unsigned int j = 0; j < gradCol.size(); j++) {
            BOOST_TEST(gradCol[j] == dGradChidGradPhi[j][i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_generate_input_variable_string) {
    // Initialize the time
    std::vector<double> time = {10., 2.5};

    // Initialize the material parameters
    std::vector<double> fparams = {
        2,     2.4e2,  1.5e1,                        // Macro hardening parameters
        2,     1.4e2,  2.0e1,                        // Micro hardening parameters
        2,     2.0e0,  2.7e1,                        // Micro gradient hardening parameters
        2,     0.56,   0.2,                          // Macro flow parameters
        2,     0.15,   -0.2,                         // Micro flow parameters
        2,     0.82,   0.1,                          // Micro gradient flow parameters
        2,     0.70,   0.3,                          // Macro yield parameters
        2,     0.40,   -0.3,                         // Micro yield parameters
        2,     0.52,   0.4,                          // Micro gradient yield parameters
        2,     696.47, 65.84,                        // A stiffness tensor parameters
        5,     -7.69,  -51.92, 38.61, -27.31, 5.13,  // B stiffness tensor parameters
        11,    1.85,   -0.19,  -1.08, -1.57,  2.29,
        -0.61, 5.97,   -2.02,  2.38,  -0.32,  -3.25,  // C stiffness tensor parameters
        2,     -51.92, 5.13,                          // D stiffness tensor parameters
        0.4,   0.3,    0.35,   1e-8,  1e-8            // Integration parameters
    };

    // Initialize the gradient of the macro displacement
    double current_grad_u[3][3] = {
        {0.200, 0.100, 0.000},
        {0.100, 0.001, 0.000},
        {0.000, 0.000, 0.000}
    };

    double previous_grad_u[3][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };
    // Initialize the micro displacement
    double current_phi[9] = {0.100, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

    double previous_phi[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Initialize the gradient of the micro displacement
    double current_grad_phi[9][3] = {
        {0.13890017,  -0.3598602,  -0.08048856},
        {-0.18572739, 0.06847269,  0.22931628 },
        {-0.01829735, -0.48731265, -0.25277529},
        {0.26626212,  0.4844646,   -0.31965177},
        {0.49197846,  0.19051656,  -0.0365349 },
        {-0.06607774, -0.33526875, -0.15803078},
        {0.09738707,  -0.49482218, -0.39584868},
        {-0.45599864, 0.08585038,  -0.09432794},
        {0.23055539,  0.07564162,  0.24051469 }
    };

    double previous_grad_phi[9][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    // Initialize the state variable vector
    std::vector<double> SDVS(55, 0);

    // Initialize the additional degree of freedom vectors
    std::vector<double>               current_ADD_DOF;
    std::vector<std::vector<double> > current_ADD_grad_DOF;

    std::vector<double>               previous_ADD_DOF;
    std::vector<std::vector<double> > previous_ADD_grad_DOF;

    std::string result;

    std::string answer =
        "time:\n 10.000000, 2.500000,\nfparams:\n 2.000000, 240.000000, 15.000000, 2.000000, 140.000000, 20.000000, "
        "2.000000, 2.000000, 27.000000, 2.000000, 0.560000, 0.200000, 2.000000, 0.150000, -0.200000, 2.000000, "
        "0.820000, 0.100000, 2.000000, 0.700000, 0.300000, 2.000000, 0.400000, -0.300000, 2.000000, 0.520000, "
        "0.400000, 2.000000, 696.470000, 65.840000, 5.000000, -7.690000, -51.920000, 38.610000, -27.310000, 5.130000, "
        "11.000000, 1.850000, -0.190000, -1.080000, -1.570000, 2.290000, -0.610000, 5.970000, -2.020000, 2.380000, "
        "-0.320000, -3.250000, 2.000000, -51.920000, 5.130000, 0.400000, 0.300000, 0.350000, 0.000000, "
        "0.000000,\ncurrent_grad_u:\n 0.200000, 0.100000, 0.000000,\n 0.100000, 0.001000, 0.000000,\n 0.000000, "
        "0.000000, 0.000000,\n\ncurrent_phi:\n 0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000,\ncurrent_grad_phi:\n 0.138900, -0.359860, -0.080489,\n -0.185727, 0.068473, 0.229316,\n "
        "-0.018297, -0.487313, -0.252775,\n 0.266262, 0.484465, -0.319652,\n 0.491978, 0.190517, -0.036535,\n "
        "-0.066078, -0.335269, -0.158031,\n 0.097387, -0.494822, -0.395849,\n -0.455999, 0.085850, -0.094328,\n "
        "0.230555, 0.075642, 0.240515,\n\nprevious_grad_u:\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, "
        "0.000000,\n 0.000000, 0.000000, 0.000000,\n\nprevious_phi:\n 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000,\nprevious_grad_phi:\n 0.000000, 0.000000, 0.000000,\n "
        "0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, "
        "0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, "
        "0.000000,\n 0.000000, 0.000000, 0.000000,\n\nSDVS:\n 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000, 0.000000, 0.000000, 0.000000, 0.000000, "
        "0.000000,\ncurrent_ADD_DOF:\n\ncurrent_ADD_grad_DOF:\n\nprevious_ADD_DOF:\n\nprevious_ADD_grad_DOF:\n";

    tardigradeMicromorphicLinearElasticity::generate_input_variable_string(
        time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi,
        SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, result);

    BOOST_CHECK(answer.compare(result) == 0);
}

BOOST_AUTO_TEST_CASE(testEvaluateHydraModel) {
    /*!
     * Test the model evaluation interface for the stresses computed
     * in the reference configuration.
     *
     */

    const std::vector<double> time = {10, 2.7};

    const std::vector<double> fparams = {2,  1.7, 1.8, 5,  2.8, .76, .15, 9.8, 5.4, 11, 1.,  2.,
                                         3., 4.,  5.,  6., 7.,  8.,  9.,  10., 11., 2,  .76, 5.4};

    const double current_grad_u[3][3] = {
        {-1.07901185, -1.09656192, -0.04629144},
        {-0.77749189, -1.27877771, -0.82648234},
        {0.66484637,  -0.05552567, -1.65125738}
    };

    const double current_phi[9] = {-1.40391532, -0.42715691, 0.75393369, 0.2849511,  -2.06484257,
                                   -0.52190902, 1.07238446,  0.19155907, -0.39704566};

    const double current_grad_phi[9][3] = {
        {0.14940184, 0.12460812, 0.31971128},
        {0.67550862, 0.61095383, 0.87972732},
        {0.30872424, 0.32158187, 0.25480281},
        {0.45570006, 0.69090695, 0.72388584},
        {0.14880964, 0.67520596, 0.15106516},
        {0.77810545, 0.07641724, 0.09367471},
        {0.15905979, 0.0651695,  0.52150417},
        {0.91873444, 0.5622355,  0.50199447},
        {0.26729942, 0.89858519, 0.09043229}
    };

    const double previous_grad_u[3][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    const double previous_phi[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    const double previous_grad_phi[9][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    std::vector<double>                     SDVS;
    const std::vector<double>               current_ADD_DOF, previous_ADD_DOF;
    const std::vector<std::vector<double> > current_ADD_grad_DOF, previous_ADD_grad_DOF;

    std::vector<double> PK2_answer = {-26.78976487, 91.99831835,  135.04096376, -63.68792655, 149.68226149,
                                      186.67587146, -42.54105342, 125.2317492,  150.55767059};

    std::vector<double> SIGMA_answer = {-47.59920949, 20.84881327, 93.02392773, 20.84881327, 302.43209139,
                                        311.0104045,  93.02392773, 311.0104045, 312.60512922};

    std::vector<double> M_answer = {-50.37283054, -23.25778149, -37.92963077, -19.16962188, -32.97279228, -14.89104497,
                                    -33.4026237,  -15.47947779, -40.31460994, -16.29637436, -36.63942799, -18.22777296,
                                    -39.33546661, -86.69472439, -59.29150146, -15.76480164, -55.42039768, -35.09720118,
                                    -28.94394503, -17.96726082, -45.09734176, -16.46568416, -50.79898863, -39.19129183,
                                    -47.46372724, -42.98201472, -45.57864883};

    std::vector<double>               PK2_result, SIGMA_result, M_result;
    std::vector<std::vector<double> > ADD_TERMS;
    std::string                       output_message;

    int errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
        time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi,
        SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result, SIGMA_result,
        M_result, ADD_TERMS, output_message);

    BOOST_CHECK(errorCode <= 0);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(PK2_result, PK2_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(SIGMA_result, SIGMA_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(M_result, M_answer));

    // Test of the Jacobians
    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();
    ADD_TERMS.clear();

    variableMatrix result_dPK2dGradU(9, variableVector(9, 0));

    variableMatrix result_dPK2dPhi(9, variableVector(9, 0));

    variableMatrix result_dPK2dGradPhi(9, variableVector(27, 0));

    variableMatrix result_dSIGMAdGradU(9, variableVector(9, 0));

    variableMatrix result_dSIGMAdPhi(9, variableVector(9, 0));

    variableMatrix result_dSIGMAdGradPhi(9, variableVector(27, 0));

    variableMatrix result_dMdGradU(27, variableVector(9, 0));

    variableMatrix result_dMdPhi(27, variableVector(9, 0));

    variableMatrix result_dMdGradPhi(27, variableVector(27, 0));

    std::vector<variableMatrix> ADD_JACOBIANS;

    errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
        time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi,
        SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result, SIGMA_result,
        M_result, result_dPK2dGradU, result_dPK2dPhi, result_dPK2dGradPhi, result_dSIGMAdGradU, result_dSIGMAdPhi,
        result_dSIGMAdGradPhi, result_dMdGradU, result_dMdPhi, result_dMdGradPhi, ADD_TERMS, ADD_JACOBIANS,
        output_message);

    BOOST_CHECK(errorCode <= 0);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(PK2_result, PK2_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(SIGMA_result, SIGMA_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(M_result, M_answer));

    variableMatrix dPK2dGradU(9, variableVector(9, 0));

    variableMatrix dPK2dPhi(9, variableVector(9, 0));

    variableMatrix dPK2dGradPhi(9, variableVector(27, 0));

    variableMatrix dSIGMAdGradU(9, variableVector(9, 0));

    variableMatrix dSIGMAdPhi(9, variableVector(9, 0));

    variableMatrix dSIGMAdGradPhi(9, variableVector(27, 0));

    variableMatrix dMdGradU(27, variableVector(9, 0));

    variableMatrix dMdPhi(27, variableVector(9, 0));

    variableMatrix dMdGradPhi(27, variableVector(27, 0));

    variableType eps = 1e-6;

    for (unsigned int i = 0; i < 9; i++) {
        variableVector delta(9, 0);

        unsigned int row = i / 3;

        unsigned int col = i % 3;

        delta[i] = eps * std::fabs(current_grad_u[row][col]) + eps;

        variableType current_grad_u_p[3][3];
        variableType current_grad_u_m[3][3];

        for (unsigned int _i = 0; _i < 3; _i++) {
            for (unsigned int _j = 0; _j < 3; _j++) {
                current_grad_u_p[_i][_j] = current_grad_u[_i][_j] + delta[3 * _i + _j];
                current_grad_u_m[_i][_j] = current_grad_u[_i][_j] - delta[3 * _i + _j];
            }
        }

        variableVector PK2_p, PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p, M_m;

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u_p, current_phi, current_grad_phi, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_p, SIGMA_p, M_p, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u_m, current_phi, current_grad_phi, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_m, SIGMA_m, M_m, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        for (unsigned int j = 0; j < PK2_p.size(); j++) {
            dPK2dGradU[j][i] = (PK2_p[j] - PK2_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < SIGMA_p.size(); j++) {
            dSIGMAdGradU[j][i] = (SIGMA_p[j] - SIGMA_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < M_p.size(); j++) {
            dMdGradU[j][i] = (M_p[j] - M_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dPK2dGradU, result_dPK2dGradU));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dSIGMAdGradU, result_dSIGMAdGradU));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dMdGradU, result_dMdGradU));

    for (unsigned int i = 0; i < 9; i++) {
        variableVector delta(9, 0);

        delta[i] = eps * std::fabs(current_phi[i]) + eps;

        variableType current_phi_p[9];
        variableType current_phi_m[9];

        for (unsigned int _i = 0; _i < 3; _i++) {
            for (unsigned int _j = 0; _j < 3; _j++) {
                current_phi_p[3 * _i + _j] = current_phi[3 * _i + _j] + delta[3 * _i + _j];
                current_phi_m[3 * _i + _j] = current_phi[3 * _i + _j] - delta[3 * _i + _j];
            }
        }

        variableVector PK2_p, PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p, M_m;

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u, current_phi_p, current_grad_phi, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_p, SIGMA_p, M_p, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u, current_phi_m, current_grad_phi, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_m, SIGMA_m, M_m, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        for (unsigned int j = 0; j < PK2_p.size(); j++) {
            dPK2dPhi[j][i] = (PK2_p[j] - PK2_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < SIGMA_p.size(); j++) {
            dSIGMAdPhi[j][i] = (SIGMA_p[j] - SIGMA_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < M_p.size(); j++) {
            dMdPhi[j][i] = (M_p[j] - M_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dPK2dPhi, result_dPK2dPhi));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dSIGMAdPhi, result_dSIGMAdPhi));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dMdPhi, result_dMdPhi));

    for (unsigned int i = 0; i < 27; i++) {
        variableVector delta(27, 0);

        unsigned int row = i / 9;

        unsigned int col = i % 9;

        delta[i] = eps * std::fabs(current_grad_phi[row][col]) + eps;

        variableType current_grad_phi_p[9][3];
        variableType current_grad_phi_m[9][3];

        for (unsigned int _i = 0; _i < 9; _i++) {
            for (unsigned int _j = 0; _j < 3; _j++) {
                current_grad_phi_p[_i][_j] = current_grad_phi[_i][_j] + delta[3 * _i + _j];
                current_grad_phi_m[_i][_j] = current_grad_phi[_i][_j] - delta[3 * _i + _j];
            }
        }

        variableVector PK2_p, PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p, M_m;

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u, current_phi, current_grad_phi_p, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_p, SIGMA_p, M_p, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
            time, fparams, current_grad_u, current_phi, current_grad_phi_m, previous_grad_u, previous_phi,
            previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
            PK2_m, SIGMA_m, M_m, ADD_TERMS, output_message);

        BOOST_CHECK(errorCode <= 0);

        for (unsigned int j = 0; j < PK2_p.size(); j++) {
            dPK2dGradPhi[j][i] = (PK2_p[j] - PK2_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < SIGMA_p.size(); j++) {
            dSIGMAdGradPhi[j][i] = (SIGMA_p[j] - SIGMA_m[j]) / (2 * delta[i]);
        }

        for (unsigned int j = 0; j < M_p.size(); j++) {
            dMdGradPhi[j][i] = (M_p[j] - M_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dPK2dGradPhi, result_dPK2dGradPhi));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dSIGMAdGradPhi, result_dSIGMAdGradPhi));
    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(dMdGradPhi, result_dMdGradPhi));
}

!> \file
!> $Id$
!> \brief This module handles all fmm routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all fmm class routines.
MODULE FMM_ROUTINES

  USE BASE_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE HAMILTON_JACOBI_EQUATIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  PUBLIC FMM_EQUATIONS_SET_CLASS_TYPE_GET,FMM_PROBLEM_CLASS_TYPE_GET
  
  PUBLIC FMM_EQUATIONS_SET_CLASS_TYPE_SET, &
    & FMM_EQUATIONS_SET_SETUP,FMM_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & FMM_PROBLEM_CLASS_TYPE_SET,FMM_PROBLEM_SETUP, &
    & FMM_PRE_SOLVE,FMM_POST_SOLVE, &
    & FMM_SPEED_FUNCTION_FIBER_DIRECTION
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a fmm equation set class.
  SUBROUTINE FMM_EQUATIONS_SET_CLASS_TYPE_GET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_TYPE !<On return, the equation type
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SUBTYPE !<On return, the equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("FMM_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_FMM_CLASS) THEN
        EQUATIONS_TYPE=EQUATIONS_SET%TYPE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Equations set is not the fmm type",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("FMM_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("FMM_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE FMM_EQUATIONS_SET_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a fmm equation set class.
  SUBROUTINE FMM_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_EQUATIONS_SET_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FMM_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_HJ_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_STANDARD_HJ_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("FMM_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("FMM_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE FMM_EQUATIONS_SET_CLASS_TYPE_SET


  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a fmm equations set class.
  SUBROUTINE FMM_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        CALL HJ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("FMM_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("FMM_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE FMM_EQUATIONS_SET_SETUP


  !
  !================================================================================================================================
  !

  !>Sets up the Laplace equation type of a classical field equations set class.
  SUBROUTINE HJ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("HJ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
        CALL HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a HJ equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("HJ_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("HJ_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("HJ_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE HJ_EQUATION_EQUATIONS_SET_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equation.
  SUBROUTINE HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS, component_idx, EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !CALL LAPLACE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "Activation time",ERR,ERROR,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & "Calculation state", &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & "Activation time",ERR,ERROR,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & "Calculation state",ERR,ERROR,*999)
              !!Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
                !Nothing to do
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
                !Nothing to do
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Create materials field
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              NUMBER_OF_COMPONENTS = 3
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials",ERR,ERROR,*999)

                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Parameters", &
                    & ERR,ERROR,*999)

                DO component_idx=1,NUMBER_OF_COMPONENTS
                  !Default the materials components to the geometric interpolation setup with node based interpolation
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO

                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_GET(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES)
                CASE(1)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Invalid number of variables. The number of variables for material field number "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))// &
                    & ", but should be 1"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                ! 1.49254 s/m along the fibers
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1.49254_DP,ERR,ERROR,*999)
                DO component_idx=2,NUMBER_OF_COMPONENTS
                  ! 4.25532 s/m transversal to the fibers
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,4.25532_DP,ERR,ERROR,*999)
                ENDDO
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_STATEITERATION,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_TIME_STEPPING,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FMM_SOLUTION_METHOD)
              !Nothing to do
            CASE DEFAULT
              LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE HJ_EQUATION_EQUATIONS_SET_STANDARD_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a fmm equation set class.
  SUBROUTINE FMM_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_EQUATIONS_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        !CALL HJ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        EQUATIONS_SET%SOLUTION_METHOD=SOLUTION_METHOD
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("FMM_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("FMM_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE FMM_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a fmm problem class.
  SUBROUTINE FMM_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem type
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("FMM_PROBLEM_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%CLASS==PROBLEM_FMM_CLASS) THEN
        PROBLEM_EQUATION_TYPE=PROBLEM%TYPE
        PROBLEM_SUBTYPE=PROBLEM%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Problem is not fmm class",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_PROBLEM_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("FMM_PROBLEM_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("FMM_PROBLEM_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE FMM_PROBLEM_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a fmm problem class.
  SUBROUTINE FMM_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem type
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_PROBLEM_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_EQUATION_TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        CALL HJ_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem equation type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_EQUATION_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_PROBLEM_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("FMM_PROBLEM_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("FMM_PROBLEM_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE FMM_PROBLEM_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a fmm problem class.
  SUBROUTINE FMM_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    PROCEDURE(SPEED_FUNCTION_INTERFACE), POINTER :: SPEED_FUNCTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%TYPE==PROBLEM_HJ_EQUATION_TYPE .AND. PROBLEM%SUBTYPE==PROBLEM_STANDARD_HJ_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a state solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_STATE_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            SPEED_FUNCTION=>FMM_SPEED_FUNCTION_FIBER_DIRECTION
            CALL SOLVER_FMM_SPEED_FUNCTION_SET(SOLVER,SPEED_FUNCTION,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATEITERATION,ERR,ERROR,*999)
            !CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_TIME_STEPPING,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " and subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR)) &
          & //" are not valid for a problem of the class FMM."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("FMM_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("FMM_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE FMM_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a fmm problem class.
  SUBROUTINE FMM_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        !CALL HJ_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_PRE_SOLVE")
    RETURN
999 CALL ERRORS("FMM_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("FMM_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE FMM_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a fmm problem class.
  SUBROUTINE FMM_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FMM_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        !CALL HJ_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fmm problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FMM_POST_SOLVE")
    RETURN
999 CALL ERRORS("FMM_POST_SOLVE",ERR,ERROR)
    CALL EXITS("FMM_POST_SOLVE")
    RETURN 1
  END SUBROUTINE FMM_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Laplace equation type .
  SUBROUTINE HJ_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("HJ_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_HJ_SUBTYPE )        
        PROBLEM%CLASS=PROBLEM_FMM_CLASS
        PROBLEM%TYPE=PROBLEM_HJ_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_HJ_SUBTYPE     
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for HJ type of the FMM problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("HJ_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("HJ_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("HJ_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE HJ_PROBLEM_SUBTYPE_SET


  !
  !================================================================================================================================
  !
  
  SUBROUTINE FMM_SPEED_FUNCTION_FIBER_DIRECTION(EQUATIONS_SET,NODE_A,NODE_B,TIME,ERR,ERROR)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: NODE_A
    INTEGER(INTG), INTENT(IN) :: NODE_B
    REAL(DP), INTENT(OUT) :: TIME
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD, MATERIALS_FIELD, GEOMETRIC_FIELD, FIBRE_FIELD
    REAL(DP) :: time_a, time_b
    REAL(DP), ALLOCATABLE :: fibre_direction_matrix_a(:,:), fibre_direction_matrix_b(:,:)
    REAL(DP), ALLOCATABLE :: fibre_direction_a(:), fibre_direction_b(:)
    REAL(DP), ALLOCATABLE :: speed_a(:), speed_b(:)
    REAL(DP), ALLOCATABLE :: dist(:)
    INTEGER(INTG) :: d
    REAL(DP), ALLOCATABLE :: coordinates_a(:), coordinates_b(:)
    INTEGER(INTG) :: num_dims
    
    CALL ENTERS("FMM_SPEED_FUNCTION_FIBER_DIRECTION",ERR,ERROR,*999)
    
    num_dims = 3
    
    ALLOCATE(coordinates_a(num_dims),STAT=ERR)
    ALLOCATE(coordinates_b(num_dims),STAT=ERR)
    ALLOCATE(dist(num_dims),STAT=ERR)
    ALLOCATE(fibre_direction_a(num_dims),STAT=ERR)
    ALLOCATE(fibre_direction_b(num_dims),STAT=ERR)
    ALLOCATE(fibre_direction_matrix_a(num_dims,num_dims),STAT=ERR)
    ALLOCATE(fibre_direction_matrix_b(num_dims,num_dims),STAT=ERR)
    ALLOCATE(speed_a(num_dims),STAT=ERR)
    ALLOCATE(speed_b(num_dims),STAT=ERR)
    
    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
    MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
    GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
    FIBRE_FIELD=>EQUATIONS_SET%GEOMETRY%FIBRE_FIELD
    
    DO d=1,num_dims
      ! Get coordinates from geometric field
      CALL FIELD_PARAMETER_SET_GET_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_A, &
        & d,coordinates_a(d),ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_B, &
        & d,coordinates_b(d),ERR,ERROR,*999)
      dist(d) = coordinates_b(d)-coordinates_a(d)
      
      ! Get fibre orientation from fibre field
      CALL FIELD_PARAMETER_SET_GET_NODE(FIBRE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_A, &
        & d,fibre_direction_a(d),ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(FIBRE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_B, &
        & d,fibre_direction_b(d),ERR,ERROR,*999)
      
      ! Get speed from material field
      CALL FIELD_PARAMETER_SET_GET_NODE(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_A, &
        & d,speed_a(d),ERR,ERROR,*999)
      CALL FIELD_PARAMETER_SET_GET_NODE(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,NODE_B, &
        & d,speed_b(d),ERR,ERROR,*999)
    ENDDO
    
    CALL FMM_CREATE_FIBER_MATRIX_TRANSPOSED(fibre_direction_a, fibre_direction_matrix_a, ERR, ERROR, *999)
    CALL FMM_CALCULATE_TRAVEL_TIME(fibre_direction_matrix_a, dist, speed_a, time_a, ERR, ERROR, *999)
    
    CALL FMM_CREATE_FIBER_MATRIX_TRANSPOSED(fibre_direction_b, fibre_direction_matrix_b, ERR, ERROR, *999)
    CALL FMM_CALCULATE_TRAVEL_TIME(fibre_direction_matrix_b, dist, speed_b, time_b, ERR, ERROR, *999)
    
    TIME = (TIME_A + TIME_B) / 2.0_DP
    
    DEALLOCATE(coordinates_a)
    DEALLOCATE(coordinates_b)
    DEALLOCATE(dist)
    DEALLOCATE(fibre_direction_a)
    DEALLOCATE(fibre_direction_b)
    DEALLOCATE(fibre_direction_matrix_a)
    DEALLOCATE(fibre_direction_matrix_b)
    DEALLOCATE(speed_a)
    DEALLOCATE(speed_b)
    
    CALL EXITS("FMM_SPEED_FUNCTION_FIBER_DIRECTION")
    RETURN
999 CALL ERRORS("FMM_SPEED_FUNCTION_FIBER_DIRECTION",ERR,ERROR)    
    CALL EXITS("FMM_SPEED_FUNCTION_FIBER_DIRECTION")
    RETURN 1
    
  END SUBROUTINE
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FMM_CALCULATE_TRAVEL_TIME(FIBER_DIRECTION_MATRIX,DIST,SPEED,TRAVEL_TIME,ERR,ERROR,*)
    !Argument variables
    REAL(DP), ALLOCATABLE :: FIBER_DIRECTION_MATRIX(:,:)
    REAL(DP), ALLOCATABLE :: DIST(:)
    REAL(DP), ALLOCATABLE :: SPEED(:)
    REAL(DP), INTENT(OUT) :: TRAVEL_TIME
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), DIMENSION(3) :: x
    REAL(DP) :: vmv
    
    CALL ENTERS("FMM_CALCULATE_TRAVEL_TIME",ERR,ERROR,*999)
    
    x(1) = DIST(1) * FIBER_DIRECTION_MATRIX(1,1) + DIST(2) * FIBER_DIRECTION_MATRIX(1,2) &
        & + DIST(3) * FIBER_DIRECTION_MATRIX(1,3)
    x(2) = DIST(1) * FIBER_DIRECTION_MATRIX(2,1) + DIST(2) * FIBER_DIRECTION_MATRIX(2,2) &
        & + DIST(3) * FIBER_DIRECTION_MATRIX(2,3)
    x(3) = DIST(1) * FIBER_DIRECTION_MATRIX(3,1) + DIST(2) * FIBER_DIRECTION_MATRIX(3,2) &
        & + DIST(3) * FIBER_DIRECTION_MATRIX(3,3)
    vmv = x(1) * x(1) * SPEED(1) * SPEED(1) + x(2) * x(2) * SPEED(2) * SPEED(2) + x(3) * x(3) * SPEED(3) * SPEED(3)
    TRAVEL_TIME = SQRT(vmv)
    
    CALL EXITS("FMM_CALCULATE_TRAVEL_TIME")
    RETURN
999 CALL ERRORS("FMM_CALCULATE_TRAVEL_TIME",ERR,ERROR)    
    CALL EXITS("FMM_CALCULATE_TRAVEL_TIME")
    RETURN 1
    
  END SUBROUTINE FMM_CALCULATE_TRAVEL_TIME
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FMM_CREATE_FIBER_MATRIX_TRANSPOSED(FIBER_DIRECTION,MATRIX,ERR,ERROR,*)
    !Argument variables
    REAL(DP), ALLOCATABLE :: FIBER_DIRECTION(:)
    REAL(DP), ALLOCATABLE :: MATRIX(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: normInv, projection, normInvSecond
    REAL(DP), DIMENSION(3) :: orthoVector

    CALL ENTERS("FMM_CREATE_FIBER_MATRIX_TRANSPOSED",ERR,ERROR,*999)
    
    ! 1st vector: normalized fiber direction
    normInv = 1.0_DP / SQRT(FIBER_DIRECTION(1) * FIBER_DIRECTION(1) + FIBER_DIRECTION(2) * FIBER_DIRECTION(2) &
        & + FIBER_DIRECTION(3) * FIBER_DIRECTION(3))
    MATRIX(1,1) = FIBER_DIRECTION(1) * normInv
    MATRIX(1,2) = FIBER_DIRECTION(2) * normInv
    MATRIX(1,3) = FIBER_DIRECTION(3) * normInv
    IF (ABS(normInv - 1.0_DP) > 0.001_DP) THEN
      !TODO: Give a warning here
    ENDIF
    
    ! 2nd and 3rd vector: create an orthonormal system
    
    ! 2nd vector: create a non-collinear vector, orthogonalize it to the first one and normalize it
    IF (ABS(FIBER_DIRECTION(1)) < 0.5) THEN
      orthoVector(1) = 1
      orthoVector(2) = 0
      orthoVector(3) = 0
    ELSE
      orthoVector(1) = 0
      orthoVector(2) = 1
      orthoVector(3) = 0
    ENDIF
    projection = orthoVector(1) * MATRIX(1,1) + orthoVector(2) * MATRIX(1,2) + orthoVector(3) * MATRIX(1,3)
    orthoVector(1) = orthoVector(1) - projection * MATRIX(1,1)
    orthoVector(2) = orthoVector(2) - projection * MATRIX(1,2)
    orthoVector(3) = orthoVector(3) - projection * MATRIX(1,3)
    normInvSecond = 1.0_DP / SQRT(orthoVector(1) * orthoVector(1) + orthoVector(2) * orthoVector(2) &
        & + orthoVector(3) * orthoVector(3))
    MATRIX(2,1) = orthoVector(1) * normInvSecond
    MATRIX(2,2) = orthoVector(2) * normInvSecond
    MATRIX(2,3) = orthoVector(3) * normInvSecond
    
    ! 3rd vector: cross product of the first two
    MATRIX(3,1) = MATRIX(1,2) * MATRIX(2,3) - MATRIX(1,3) * MATRIX(2,2)
    MATRIX(3,2) = MATRIX(1,3) * MATRIX(2,1) - MATRIX(1,1) * MATRIX(2,3)
    MATRIX(3,3) = MATRIX(1,1) * MATRIX(2,2) - MATRIX(1,2) * MATRIX(2,1)
    
    CALL EXITS("FMM_CREATE_FIBER_MATRIX_TRANSPOSED")
    RETURN
999 CALL ERRORS("FMM_CREATE_FIBER_MATRIX_TRANSPOSED",ERR,ERROR)    
    CALL EXITS("FMM_CREATE_FIBER_MATRIX_TRANSPOSED")
    RETURN 1
    
  END SUBROUTINE FMM_CREATE_FIBER_MATRIX_TRANSPOSED
  
  !
  !================================================================================================================================
  !
  
END MODULE FMM_ROUTINES


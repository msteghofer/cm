!> \file
!> $Id$
!> \author Ali Pashaei and Martin Steghöfer
!> \brief This module handles all Mesh I/O routines.
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

!>This module handles all Hamilton-Jacobi equations routines.
MODULE MESH_IO_ROUTINES
  
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES
  
  IMPLICIT NONE
  
  PRIVATE
  
  !Module parameters
  
  !Module types
  TYPE VTK_FIELD_TYPE
    REAL(DP), ALLOCATABLE :: VALUES(:)
  END TYPE VTK_FIELD_TYPE
  
  TYPE VTK_FIELDS_SET_TYPE
    TYPE(VTK_FIELD_TYPE), ALLOCATABLE :: POINT_ASSOCIATED_FIELDS(:)
    TYPE(VTK_FIELD_TYPE), ALLOCATABLE :: CELL_ASSOCIATED_FIELDS(:)
  END TYPE VTK_FIELDS_SET_TYPE
  
  !Module variables
  
  !Interfaces
  
  PUBLIC READ_TETGEN_MESH,READ_VTK_MESH,WRITE_VTK_MESH
  
CONTAINS

#define MAX_CHARS_PER_LINE 150
#define VTK_EMPTY_CELL       0
#define VTK_VERTEX           1
#define VTK_POLY_VERTEX      2
#define VTK_LINE             3
#define VTK_POLY_LINE        4
#define VTK_TRIANGLE         5
#define VTK_TRIANGLE_STRIP   6
#define VTK_POLYGON          7
#define VTK_PIXEL            8
#define VTK_QUAD             9
#define VTK_TETRA            10
#define VTK_VOXEL            11
#define VTK_HEXAHEDRON       12
#define VTK_WEDGE            13
#define VTK_PYRAMID          14
#define VTK_PENTAGONAL_PRISM 15
#define VTK_HEXAGONAL_PRISM  16

  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_TETGEN_MESH(INPUT_FILE_NAME,WORLD_REGION,REGION_NUMBER,MESH,GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VARYING_STRING), INTENT(IN) :: INPUT_FILE_NAME
    TYPE(REGION_TYPE), INTENT(IN), POINTER :: WORLD_REGION
    INTEGER(INTG), INTENT(IN) :: REGION_NUMBER
    TYPE(MESH_TYPE), INTENT(INOUT), POINTER :: MESH
    TYPE(FIELD_TYPE), INTENT(INOUT), POINTER :: GEOMETRIC_FIELD
    TYPE(DECOMPOSITION_TYPE), INTENT(INOUT), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    INTEGER(INTG) :: TEXT_LENGTH, NUMBER_OF_NODES, NUMBER_OF_ELEMENTS, NUMBER_OF_NODE_COMPONENTS, &
        & NUMBER_OF_ELEMENT_COMPONENTS, ELEMENT_ID, NODE_ID, NUMBER_OF_PROCESSORS, TEMP_INT, i, j
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    REAL(DP), ALLOCATABLE :: NODE_COORDINATES(:)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("READ_TETGEN_MESH",ERR,ERROR,*999)
    
    OPEN (11,FILE=CHAR(INPUT_FILE_NAME//".node"))
    READ (11,*) NUMBER_OF_NODES, NUMBER_OF_NODE_COMPONENTS, TEMP_INT, TEMP_INT
    
    OPEN (12,FILE=CHAR(INPUT_FILE_NAME//".ele"))
    READ (12,*) NUMBER_OF_ELEMENTS, NUMBER_OF_ELEMENT_COMPONENTS, TEMP_INT
    
    ! Create coordinate system
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_CREATE_START(REGION_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)
    
    ! Create region and assign coordinate system to it
    NULLIFY(REGION)
    CALL REGION_CREATE_START(REGION_NUMBER,WORLD_REGION,REGION,ERR,ERROR,*999)
    CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)
    
    ! Create linear basis
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START(REGION_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_SET(BASIS,BASIS_SIMPLEX_TYPE,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
        & BASIS_LINEAR_SIMPLEX_INTERPOLATION/),ERR,ERROR,*999)
    CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
    
    !Create mesh
    NULLIFY(MESH)
    CALL MESH_CREATE_START(REGION_NUMBER,REGION,NUMBER_OF_NODE_COMPONENTS,MESH,ERR,ERROR,*999)
    
    ! Create nodes
    NULLIFY(NODES)
    CALL NODES_CREATE_START(REGION,NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
    CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
    
    ! Create elements
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,1,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
    NULLIFY(ELEMENTS)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,ELEMENTS,ERR,ERROR,*999)
    ALLOCATE(ELEMENT_NODES(NUMBER_OF_ELEMENT_COMPONENTS),STAT=ERR)
    DO i=1,NUMBER_OF_ELEMENTS
      READ (12,*) ELEMENT_ID,ELEMENT_NODES
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ELEMENT_ID,ELEMENTS,ELEMENT_NODES,Err,ERROR,*999)
    ENDDO
    DEALLOCATE(ELEMENT_NODES)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*999)
    
    CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999)
    
    !Calculate decomposition
    NULLIFY(DECOMPOSITION)
    NUMBER_OF_PROCESSORS = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    CALL DECOMPOSITION_CREATE_START(REGION_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_PROCESSORS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*999)

    !Create a field to put the geometry
    NULLIFY(GEOMETRIC_FIELD)
    CALL FIELD_CREATE_START(REGION_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    CALL FIELD_TYPE_SET(GEOMETRIC_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_VARIABLES_SET(GEOMETRIC_FIELD,1,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_COMPONENTS_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_NODE_COMPONENTS,Err,ERROR,*999)
    DO i=1,NUMBER_OF_NODE_COMPONENTS
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,i,1,ERR,ERROR,*999)
    ENDDO
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)

    !Set node positions
    ALLOCATE(NODE_COORDINATES(NUMBER_OF_NODE_COMPONENTS),STAT=ERR)
    DO i=1,NUMBER_OF_NODES
      READ (11,*) NODE_ID,NODE_COORDINATES
      DO j=1,NUMBER_OF_NODE_COMPONENTS
        CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,i,j, &
          & NODE_COORDINATES(j),ERR,ERROR,*999)
      ENDDO
    ENDDO
    DEALLOCATE(NODE_COORDINATES)
    
    CLOSE(11)
    CLOSE(12)
    
    CALL EXITS("READ_TETGEN_MESH")
    RETURN
999 CALL ERRORS("READ_TETGEN_MESH",ERR,ERROR)
    CALL EXITS("READ_TETGEN_MESH")
    RETURN 1
    
  END SUBROUTINE READ_TETGEN_MESH
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_INTTYPE_VALUES(FILE_HANDLE,MIN_BOUND,MAX_BOUND,CHECK_BOUNDS,NUM_VALUES_TO_READ, &
        & TEMP_STRING_CLEAN,VALUES,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    INTEGER(LINTG), INTENT(IN) :: MIN_BOUND
    INTEGER(LINTG), INTENT(IN) :: MAX_BOUND
    LOGICAL, INTENT(IN) :: CHECK_BOUNDS
    INTEGER(INTG), INTENT(IN) :: NUM_VALUES_TO_READ
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    INTEGER(LINTG), ALLOCATABLE, INTENT(OUT) :: VALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING
    INTEGER(INTG) :: i, j, NUM_VALUES_READ
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_INTTYPE_VALUES",ERR,ERROR,*999)
    
    ALLOCATE(VALUES(NUM_VALUES_TO_READ))
    NUM_VALUES_READ=0
    DO WHILE (NUM_VALUES_READ < NUM_VALUES_TO_READ)
      TEMP_STRING=TEMP_STRING_CLEAN
      READ (FILE_HANDLE,'(A)') TEMP_STRING
      READ(TEMP_STRING, *, end=881, err=882) (VALUES(i), i = NUM_VALUES_READ+1, NUM_VALUES_TO_READ)
      !LOCAL_ERROR="Found more than "//TRIM(NUMBER_TO_VSTRING(NUM_VALUES_TO_READ,"*",ERR,ERROR))//" values."
      !CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881   CONTINUE
      DO j=NUM_VALUES_READ+1,i-1
        IF (CHECK_BOUNDS) THEN
          IF (VALUES(j) < MIN_BOUND) THEN
            LOCAL_ERROR="Found value "//TRIM(NUMBER_TO_VSTRING(VALUES(j),"*",ERR,ERROR))// &
                & " less than the minimum value "//TRIM(NUMBER_TO_VSTRING(MIN_BOUND,"*",ERR,ERROR))// &
                & " allowed for the specified data type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSEIF (VALUES(j) > MAX_BOUND) THEN
            LOCAL_ERROR="Found value "// &
                & TRIM(NUMBER_TO_VSTRING(VALUES(j),"*",ERR,ERROR))// &
                " greater than the maximum value "//TRIM(NUMBER_TO_VSTRING(MAX_BOUND,"*",ERR,ERROR))// &
                & " allowed for the specified data type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDDO
      NUM_VALUES_READ=i-1
    END DO
    
    CALL EXITS("READ_VTK_INTTYPE_VALUES")
    RETURN
882 LOCAL_ERROR='Error parsing values of line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
999 CALL ERRORS("READ_VTK_INTTYPE_VALUES",ERR,ERROR)
    CALL EXITS("READ_VTK_INTTYPE_VALUES")
    RETURN 1
    
  END SUBROUTINE READ_VTK_INTTYPE_VALUES
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_FLOATTYPE_VALUES(FILE_HANDLE,MIN_BOUND,MAX_BOUND,CHECK_BOUNDS,NUM_VALUES_TO_READ, &
        & TEMP_STRING_CLEAN,VALUES,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    REAL(DP), INTENT(IN) :: MIN_BOUND
    REAL(DP), INTENT(IN) :: MAX_BOUND
    LOGICAL, INTENT(IN) :: CHECK_BOUNDS
    INTEGER(INTG), INTENT(IN) :: NUM_VALUES_TO_READ
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: VALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING
    INTEGER(INTG) :: i, j, NUM_VALUES_READ
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_FLOATTYPE_VALUES",ERR,ERROR,*999)
    
    ALLOCATE(VALUES(NUM_VALUES_TO_READ))
    NUM_VALUES_READ=0
    DO WHILE (NUM_VALUES_READ < NUM_VALUES_TO_READ)
      TEMP_STRING=TEMP_STRING_CLEAN
      READ (FILE_HANDLE,'(A)') TEMP_STRING
      READ(TEMP_STRING, *, end=881, err=882) (VALUES(i), i = NUM_VALUES_READ+1, NUM_VALUES_TO_READ)
      !LOCAL_ERROR="Found more than "//TRIM(NUMBER_TO_VSTRING(NUM_VALUES_TO_READ,"*",ERR,ERROR))//" values."
      !CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881   CONTINUE
      DO j=NUM_VALUES_READ+1,i-1
        IF (CHECK_BOUNDS) THEN
          IF (VALUES(j) < MIN_BOUND) THEN
            LOCAL_ERROR="Found value "//TRIM(NUMBER_TO_VSTRING(VALUES(j),"*",ERR,ERROR))// &
                & " less than the minimum value "//TRIM(NUMBER_TO_VSTRING(MIN_BOUND,"*",ERR,ERROR))// &
                & " allowed for the specified data type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSEIF (VALUES(j) > MAX_BOUND) THEN
            LOCAL_ERROR="Found value "// &
                & TRIM(NUMBER_TO_VSTRING(VALUES(j),"*",ERR,ERROR))// &
                " greater than the maximum value "//TRIM(NUMBER_TO_VSTRING(MAX_BOUND,"*",ERR,ERROR))// &
                & " allowed for the specified data type."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDDO
      NUM_VALUES_READ=i-1
    END DO
    
    CALL EXITS("READ_VTK_FLOATTYPE_VALUES")
    RETURN
882 LOCAL_ERROR='Error parsing values of line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
999 CALL ERRORS("READ_VTK_FLOATTYPE_VALUES",ERR,ERROR)
    CALL EXITS("READ_VTK_FLOATTYPE_VALUES")
    RETURN 1
    
  END SUBROUTINE READ_VTK_FLOATTYPE_VALUES
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_VALUES(FILE_HANDLE,TYPE_STRING,NUM_VALUES_TO_READ,TEMP_STRING_CLEAN,VALUES,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    TYPE(VARYING_STRING), INTENT(IN) :: TYPE_STRING
    INTEGER(INTG), INTENT(IN) :: NUM_VALUES_TO_READ
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: VALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    INTEGER(INTG) :: NUMBER_OF_NODES
    INTEGER(LINTG), ALLOCATABLE :: VALUES_INT(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_VALUES",ERR,ERROR,*999)
    
    SELECT CASE(CHAR(TYPE_STRING))
      CASE("bit")
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,1_LINTG,.True.,NUM_VALUES_TO_READ, &
            & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("unsigned_char")
        LOCAL_ERROR='Data type "unsigned char" not supported by OpenCMISS.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      CASE("char")
        LOCAL_ERROR='Data type "char" not supported by OpenCMISS.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      CASE("unsigned_short")
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,65535_LINTG,.True.,NUM_VALUES_TO_READ, &
            & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("short")
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,-32768_LINTG,32767_LINTG,.True.,NUM_VALUES_TO_READ, &
            & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("unsigned_int")
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,4294967295_LINTG,.True.,NUM_VALUES_TO_READ, &
            & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("int")
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,-2147483648_LINTG,2147483647_LINTG,.True.,NUM_VALUES_TO_READ, &
            & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("unsigned_long")
        !CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,18446744073709551615_LINTG,.True.,NUM_VALUES_TO_READ, &
        !    & TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        LOCAL_ERROR='Data type "unsigned long" not supported by OpenCMISS.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      CASE("long")
        !TODO: Lower bound should be -9223372036854775808, but that isn't accepted by the compiler (why? it *should* fit in an 8 byte integer!)
        CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,-9223372036854775807_LINTG,9223372036854775807_LINTG,.True., &
            & NUM_VALUES_TO_READ,TEMP_STRING_CLEAN,VALUES_INT,ERR,ERROR,*999)
        ALLOCATE(VALUES(NUM_VALUES_TO_READ))
        VALUES=VALUES_INT
        DEALLOCATE(VALUES_INT)
      CASE("float")
        CALL READ_VTK_FLOATTYPE_VALUES(FILE_HANDLE,-3.4028234663852886E38_DP,3.4028234663852886E38_DP,.True., &
            & NUM_VALUES_TO_READ,TEMP_STRING_CLEAN,VALUES,ERR,ERROR,*999)
      CASE("double")
        CALL READ_VTK_FLOATTYPE_VALUES(FILE_HANDLE,0.0_DP,0.0_DP,.False., &
            & NUM_VALUES_TO_READ,TEMP_STRING_CLEAN,VALUES,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR='Unknown data type "'//TYPE_STRING//'"!'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    
    CALL EXITS("READ_VTK_VALUES")
    RETURN
999 CALL ERRORS("READ_VTK_VALUES",ERR,ERROR)
    CALL EXITS("READ_VTK_VALUES")
    RETURN 1
    
  END SUBROUTINE READ_VTK_VALUES
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_POINTS(FILE_HANDLE,FIRST_LINE,TEMP_STRING_CLEAN,NODE_COORDINATES,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: FIRST_LINE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: NODE_COORDINATES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    INTEGER(INTG) :: NUMBER_OF_NODES
    TYPE(VARYING_STRING) :: TYPE_STRING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_POINTS",ERR,ERROR,*999)
    
    ! Check, if there are more then 2 parameters
    TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_3=TEMP_STRING_CLEAN
    READ (FIRST_LINE,*,end=881) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    ! Still here (no jump to 881) => TEMP_STRING_WORD_4 (3rd parameter) was found => Error
    LOCAL_ERROR="POINTS section keyword only takes 2 parameters, "// &
        & 'but found 3rd parameter "'//TRIM(TEMP_STRING_WORD_4)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881 CONTINUE ! everything OK, TEMP_STRING_WORD_4 (3rd parameter) was not found
    ! Read the two parameters
    READ (FIRST_LINE,*,end=882) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
    IF (.FALSE.) THEN
882   LOCAL_ERROR='POINTS keyword expects 2 parameters, but less than 2 given in line "'// &
          & TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    ! Parse integer NUMBER_OF_NODES
    READ (TEMP_STRING_WORD_2,*,end=883,err=883) NUMBER_OF_NODES
    IF (.FALSE.) THEN
883   LOCAL_ERROR='Unable to read the number of nodes in POINTS section parameters, maybe the given number "'// &
          & TRIM(TEMP_STRING_WORD_2)//'" is not an integer expresseion.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    TYPE_STRING=TEMP_STRING_WORD_3
    CALL READ_VTK_VALUES(FILE_HANDLE,TYPE_STRING,NUMBER_OF_NODES*3,TEMP_STRING_CLEAN, &
        & NODE_COORDINATES,ERR,ERROR,*999)
    
    CALL EXITS("READ_VTK_POINTS")
    RETURN
999 CALL ERRORS("READ_VTK_POINTS",ERR,ERROR)
    CALL EXITS("READ_VTK_POINTS")
    RETURN 1
    
  END SUBROUTINE READ_VTK_POINTS
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_CELLS(FILE_HANDLE,FIRST_LINE,TEMP_STRING_CLEAN,CELLS,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: FIRST_LINE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    INTEGER(LINTG), ALLOCATABLE, INTENT(OUT) :: CELLS(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    INTEGER(INTG) :: NUMBER_OF_CELLS, NUMBER_OF_CELL_VALUES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_CELLS",ERR,ERROR,*999)
    
    ! Check, if there are more then 2 parameters
    TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_3=TEMP_STRING_CLEAN
    READ (FIRST_LINE,*,end=881) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    !Still here (no jump to 881) => TEMP_STRING_WORD_4 (3rd parameter) was found => Error
    LOCAL_ERROR="CELLS section keyword only takes 2 parameters, "// &
        & 'but found 3rd parameter "'//TRIM(TEMP_STRING_WORD_4)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881 CONTINUE !everything OK, TEMP_STRING_WORD_4 (3rd parameter) was not found
    ! Read the two parameters
    READ (FIRST_LINE,*,end=882) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
    IF (.FALSE.) THEN
882   LOCAL_ERROR='CELLS keyword expects 2 parameters, but less than 2 given in line "'// &
          & TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    ! Parse integer NUMBER_OF_CELLS
    READ (TEMP_STRING_WORD_2,*,end=883,err=883) NUMBER_OF_CELLS
    IF (.FALSE.) THEN
883   LOCAL_ERROR='Unable to read the number of cells in CELLS section parameters, maybe the given number "'// &
          & TRIM(TEMP_STRING_WORD_2)//'" is not an integer expresseion.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    ! Parse integer NUMBER_OF_CELL_VALUES
    READ (TEMP_STRING_WORD_3,*,end=884,err=884) NUMBER_OF_CELL_VALUES
    IF (.FALSE.) THEN
884   LOCAL_ERROR='Unable to read the number of values in CELLS section parameters, maybe the given number "'// &
          & TRIM(TEMP_STRING_WORD_3)//'" is not an integer expresseion.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,0_LINTG,.False.,NUMBER_OF_CELL_VALUES, &
            & TEMP_STRING_CLEAN,CELLS,ERR,ERROR,*999)
    
    CALL EXITS("READ_VTK_CELLS")
    RETURN
999 CALL ERRORS("READ_VTK_CELLS",ERR,ERROR)
    CALL EXITS("READ_VTK_CELLS")
    RETURN 1
    
  END SUBROUTINE READ_VTK_CELLS
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_CELL_TYPES(FILE_HANDLE,FIRST_LINE,TEMP_STRING_CLEAN,CELL_TYPES,ERR,ERROR,*)
    ! subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: FIRST_LINE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    INTEGER(LINTG), ALLOCATABLE, INTENT(OUT) :: CELL_TYPES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
    INTEGER(INTG) :: NUMBER_OF_CELLS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_CELL_TYPES",ERR,ERROR,*999)
    
    ! Check, if there is more then 1 parameter
    TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
    READ (FIRST_LINE,*,end=881) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
    !Still here (no jump to 881) => TEMP_STRING_WORD_3 (2nd parameter) was found => Error
    LOCAL_ERROR="CELL_TYPES section keyword only takes 1 parameter, "// &
        & 'but found 2nd parameter "'//TRIM(TEMP_STRING_WORD_3)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881 CONTINUE !everything OK, TEMP_STRING_WORD_3 (2nd parameter) was not found
    ! Read the one parameter
    READ (FIRST_LINE,*,end=882) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2
    IF (.FALSE.) THEN
882   LOCAL_ERROR='CELL_TYPES keyword expects 1 parameter, but less than 1 given in line "'// &
          & TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    ! Parse integer NUMBER_OF_CELLS
    READ (TEMP_STRING_WORD_2,*,end=883,err=883) NUMBER_OF_CELLS
    IF (.FALSE.) THEN
883   LOCAL_ERROR='Unable to read the number of cells in CELL_TYPES section parameters, '// &
          & 'maybe the given number "'//TRIM(TEMP_STRING_WORD_2)//'" is not an integer expresseion.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    CALL READ_VTK_INTTYPE_VALUES(FILE_HANDLE,0_LINTG,0_LINTG,.False.,NUMBER_OF_CELLS, &
            & TEMP_STRING_CLEAN,CELL_TYPES,ERR,ERROR,*999)
    
    CALL EXITS("READ_VTK_CELL_TYPES")
    RETURN
999 CALL ERRORS("READ_VTK_CELL_TYPES",ERR,ERROR)
    CALL EXITS("READ_VTK_CELL_TYPES")
    RETURN 1
    
  END SUBROUTINE READ_VTK_CELL_TYPES
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_INCREMENT_NUMBER_OF_FIELDS(FIELD_DATA,IN_CELL_DATA,NUMBER_OF_FIELDS,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VTK_FIELDS_SET_TYPE), INTENT(INOUT) :: FIELD_DATA
    LOGICAL, INTENT(IN) :: IN_CELL_DATA
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_FIELDS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(VTK_FIELD_TYPE), ALLOCATABLE :: NEW_FIELDS(:), OLD_FIELDS(:)
    INTEGER(INTG) :: i
    
    IF (NUMBER_OF_FIELDS < 0) THEN
      LOCAL_ERROR='Cannot reduce size of fields array.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    IF (NUMBER_OF_FIELDS > 0) THEN
      ! get old fields array
      IF (IN_CELL_DATA) THEN
        CALL MOVE_ALLOC(FIELD_DATA%CELL_ASSOCIATED_FIELDS,OLD_FIELDS)
      ELSE
        CALL MOVE_ALLOC(FIELD_DATA%POINT_ASSOCIATED_FIELDS,OLD_FIELDS)
      ENDIF
      ! create new (longer) fields array
      ALLOCATE(NEW_FIELDS(SIZE(OLD_FIELDS,1)+NUMBER_OF_FIELDS))
      ! copy old fields into new array
      DO i=1,SIZE(OLD_FIELDS,1)
        NEW_FIELDS(i)=OLD_FIELDS(i)
      ENDDO
      ! write new array back to initial data structure
      IF (IN_CELL_DATA) THEN
        CALL MOVE_ALLOC(NEW_FIELDS,FIELD_DATA%CELL_ASSOCIATED_FIELDS)
      ELSE
        CALL MOVE_ALLOC(NEW_FIELDS,FIELD_DATA%POINT_ASSOCIATED_FIELDS)
      ENDIF
    ENDIF
    
    CALL EXITS("READ_VTK_INCREMENT_NUMBER_OF_FIELDS")
    RETURN
999 CALL ERRORS("READ_VTK_INCREMENT_NUMBER_OF_FIELDS",ERR,ERROR)
    CALL EXITS("READ_VTK_INCREMENT_NUMBER_OF_FIELDS")
    RETURN 1
    
  END SUBROUTINE READ_VTK_INCREMENT_NUMBER_OF_FIELDS
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_FIELD(FILE_HANDLE,FIRST_LINE,TEMP_STRING_CLEAN,FIELD_DATA,IN_CELL_DATA,ERR,ERROR,*)
    !subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: FIRST_LINE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    TYPE(VTK_FIELDS_SET_TYPE), INTENT(INOUT) :: FIELD_DATA
    LOGICAL, INTENT(IN) :: IN_CELL_DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, &
        & TEMP_STRING_WORD_4, TEMP_STRING_WORD_5
    LOGICAL :: EOF
    INTEGER(INTG) :: NUMBER_OF_FIELDS, i, j, TUPLE_DIMENSION, NUMBER_OF_TUPLES
    TYPE(VARYING_STRING) :: LOCAL_ERROR, TYPE_NAME
    
    CALL ENTERS("READ_VTK_FIELD",ERR,ERROR,*999)
    
    TEMP_STRING=FIRST_LINE
    READ (TEMP_STRING,*,end=881) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    !Still here (no jump to 881) => TEMP_STRING_WORD_4 was found => Error
    LOCAL_ERROR='FIELD only accepts 2 parameters, but extra (3rd) parameter "'// &
        & TRIM(TEMP_STRING_WORD_4)//'" was found in line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881 CONTINUE !everything OK, TEMP_STRING_WORD_4 was not found
    TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_3=TEMP_STRING_CLEAN
    READ (TEMP_STRING,*,end=882) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
    READ (TEMP_STRING_WORD_3,*,err=883) NUMBER_OF_FIELDS
    IF (NUMBER_OF_FIELDS < 0) THEN
      LOCAL_ERROR='Negative number of field variables ('//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_FIELDS,"*",ERR,ERROR))// &
          & ') found in line"'//TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    CALL READ_VTK_INCREMENT_NUMBER_OF_FIELDS(FIELD_DATA,IN_CELL_DATA,NUMBER_OF_FIELDS,ERR,ERROR,*999)
    DO i=1,NUMBER_OF_FIELDS
888   TEMP_STRING=TEMP_STRING_CLEAN
      READ (FILE_HANDLE,'(A)') TEMP_STRING
      TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
      READ (TEMP_STRING,*,end=888) TEMP_STRING_WORD_1 ! if line is empty, ignore it and read next line
      READ (TEMP_STRING,*,end=884) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, &
          & TEMP_STRING_WORD_4, TEMP_STRING_WORD_5
      !Still here (no jump to 884) => TEMP_STRING_WORD_5 was found => Error
      LOCAL_ERROR='Array only accepts 4 parameters, but extra (5th) parameter "'// &
          & TRIM(TEMP_STRING_WORD_5)//'" was found in line "'//TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
884   CONTINUE !everything OK, TEMP_STRING_WORD_5 was not found
      TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
      TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
      TEMP_STRING_WORD_3=TEMP_STRING_CLEAN
      TEMP_STRING_WORD_4=TEMP_STRING_CLEAN
      READ (TEMP_STRING,*,end=885) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
      READ (TEMP_STRING_WORD_2,*,err=886) TUPLE_DIMENSION
      READ (TEMP_STRING_WORD_3,*,err=887) NUMBER_OF_TUPLES
      TYPE_NAME=TRIM(TEMP_STRING_WORD_4)
      IF (IN_CELL_DATA) THEN
        j=SIZE(FIELD_DATA%CELL_ASSOCIATED_FIELDS,1)-(NUMBER_OF_FIELDS-i)
        CALL READ_VTK_VALUES(FILE_HANDLE,TYPE_NAME,NUMBER_OF_TUPLES*TUPLE_DIMENSION,TEMP_STRING_CLEAN, &
            & FIELD_DATA%CELL_ASSOCIATED_FIELDS(j)%VALUES,ERR,ERROR,*999)
      ELSE
        j=SIZE(FIELD_DATA%POINT_ASSOCIATED_FIELDS,1)-(NUMBER_OF_FIELDS-i)
        CALL READ_VTK_VALUES(FILE_HANDLE,TYPE_NAME,NUMBER_OF_TUPLES*TUPLE_DIMENSION,TEMP_STRING_CLEAN, &
            & FIELD_DATA%POINT_ASSOCIATED_FIELDS(j)%VALUES,ERR,ERROR,*999)
      ENDIF
    ENDDO
    
    CALL EXITS("READ_VTK_FIELD")
    RETURN
882 LOCAL_ERROR='FIELD needs 2 parameters, less than 2 were found in line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    RETURN
883 LOCAL_ERROR='Could not parse the number of fields "'//TRIM(TEMP_STRING_WORD_2)// &
        & '" (first parameter of FIELD) in line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    RETURN
885 LOCAL_ERROR='Array needs 4 parameters (name, tuple components, number of tuples and type), but '// &
        & 'less than 4 were found in line "'//TRIM(TEMP_STRING)//'".'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    RETURN
886 LOCAL_ERROR='Could not parse the tuple dimension "'//TRIM(TEMP_STRING_WORD_2)// &
        & '" (second token in line "'//TRIM(TEMP_STRING)//'").'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    RETURN
887 LOCAL_ERROR='Could not parse the number of tuples "'//TRIM(TEMP_STRING_WORD_3)// &
        & '" (third token in line "'//TRIM(TEMP_STRING)//'").'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("READ_VTK_FIELD",ERR,ERROR)
    CALL EXITS("READ_VTK_FIELD")
    RETURN 1
    
  END SUBROUTINE READ_VTK_FIELD
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_CELL_AND_POINT_DATA(FILE_HANDLE,FIRST_LINE,TEMP_STRING_CLEAN,FIELD_DATA,ERR,ERROR,*)
    !subroutine parameters
    INTEGER(INTG), INTENT(IN) :: FILE_HANDLE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: FIRST_LINE
    character(LEN=MAX_CHARS_PER_LINE), INTENT(IN) :: TEMP_STRING_CLEAN
    TYPE(VTK_FIELDS_SET_TYPE), INTENT(INOUT) :: FIELD_DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    LOGICAL :: IN_CELL_DATA
    TYPE(VARYING_STRING) :: LOCAL_ERROR, SECTION_NAME
    
    CALL ENTERS("READ_VTK_CELL_AND_POINT_DATA",ERR,ERROR,*999)
    
    TEMP_STRING=FIRST_LINE
    DO WHILE (.TRUE.)
      ! Check first keyword in line
      TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
      READ (TEMP_STRING,*,end=882) TEMP_STRING_WORD_1 ! If line was empty, ignore it and read next line (882)
      SELECT CASE(TRIM(TEMP_STRING_WORD_1))
        CASE("CELL_DATA")
          IN_CELL_DATA=.TRUE.
        CASE("POINT_DATA")
          IN_CELL_DATA=.FALSE.
        CASE("SCALARS")
          LOCAL_ERROR='The data type SCALARS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("VECTORS")
          LOCAL_ERROR='The data type VECTORS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("TENSORS")
          LOCAL_ERROR='The data type TENSORS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("NORMALS")
          LOCAL_ERROR='The data type NORMALS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("TEXTURE_COORDINATES")
          LOCAL_ERROR='The data type TEXTURE_COORDINATES is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("GLOBAL_IDS")
          LOCAL_ERROR='The data type GLOBAL_IDS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("PEDIGREE_IDS")
          LOCAL_ERROR='The data type PEDIGREE_IDS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("COLOR_SCALARS")
          LOCAL_ERROR='The data type COLOR_SCALARS is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("LOOKUP_TABLE")
          LOCAL_ERROR='The data type LOOKUP_TABLE is not supported by OpenCMISS!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        CASE("FIELD")
          CALL READ_VTK_FIELD(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,FIELD_DATA,IN_CELL_DATA,ERR,ERROR,*999)
        CASE DEFAULT
          IF (IN_CELL_DATA) THEN
            SECTION_NAME="CELL_DATA"
          ELSE
            SECTION_NAME="POINT_DATA"
          ENDIF
          LOCAL_ERROR='Found unexpected keyword "'//TRIM(TEMP_STRING_WORD_1)//'" in line "'// &
              & TRIM(TEMP_STRING)//'" in '//SECTION_NAME//' section.'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
882   CONTINUE
      ! Read next line
      TEMP_STRING=TEMP_STRING_CLEAN
      READ (FILE_HANDLE,'(A)',END=883) TEMP_STRING
    ENDDO

883 CONTINUE ! Reached EOF, no more sections    
    CALL EXITS("READ_VTK_CELL_AND_POINT_DATA")
    RETURN
999 CALL ERRORS("READ_VTK_CELL_AND_POINT_DATA",ERR,ERROR)
    CALL EXITS("READ_VTK_CELL_AND_POINT_DATA")
    RETURN 1
    
  END SUBROUTINE READ_VTK_CELL_AND_POINT_DATA
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_MESH_ERROR_CELL_TYPE(CELL_TYPE_NAME,CELL_TYPE_NUMBER,CELL_ID,ERR,ERROR,*)
    TYPE(VARYING_STRING), INTENT(IN) :: CELL_TYPE_NAME !<The error string
    INTEGER(LINTG), INTENT(IN) :: CELL_TYPE_NUMBER
    INTEGER(INTG), INTENT(IN) :: CELL_ID
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_MESH_ERROR_CELL_TYPE",ERR,ERROR,*999)
    
    LOCAL_ERROR='Cell type '//TRIM(CELL_TYPE_NAME)//' ('//TRIM(NUMBER_TO_VSTRING(CELL_TYPE_NUMBER,"*",ERR,ERROR))// &
        & ') not supported by OpenCMISS (found in cell '//TRIM(NUMBER_TO_VSTRING(CELL_ID,"*",ERR,ERROR))//')!'
    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    
    CALL EXITS("READ_VTK_MESH_ERROR_CELL_TYPE")
    RETURN
999 CALL ERRORS("READ_VTK_MESH_ERROR_CELL_TYPE",ERR,ERROR)
    CALL EXITS("READ_VTK_MESH_ERROR_CELL_TYPE")
    RETURN 1
  END SUBROUTINE READ_VTK_MESH_ERROR_CELL_TYPE
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_MESH_CHECK_CELL_NODE_NUMBER(CELL_TYPE_NAME,NUMBER_OF_NODES_DESIRED,CELL_NUMBER_OF_NODES,CELL_ID,ERR,ERROR,*)
    TYPE(VARYING_STRING), INTENT(IN) :: CELL_TYPE_NAME !<The error string
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES_DESIRED
    INTEGER(LINTG), INTENT(IN) :: CELL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(IN) :: CELL_ID
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_MESH_CHECK_CELL_NODE_NUMBER",ERR,ERROR,*999)
    
    IF (NUMBER_OF_NODES_DESIRED /= CELL_NUMBER_OF_NODES) THEN
      LOCAL_ERROR='Cells of type '//TRIM(CELL_TYPE_NAME)//' must have exactly '// &
          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_NODES_DESIRED,"*",ERR,ERROR))//' nodes, but cell '// &
          & TRIM(NUMBER_TO_VSTRING(CELL_ID,"*",ERR,ERROR))//' has '// &
          & TRIM(NUMBER_TO_VSTRING(CELL_NUMBER_OF_NODES,"*",ERR,ERROR))//' nodes!'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("READ_VTK_MESH_CHECK_CELL_NODE_NUMBER")
    RETURN
999 CALL ERRORS("READ_VTK_MESH_CHECK_CELL_NODE_NUMBER",ERR,ERROR)
    CALL EXITS("READ_VTK_MESH_CHECK_CELL_NODE_NUMBER")
    RETURN 1
  END SUBROUTINE READ_VTK_MESH_CHECK_CELL_NODE_NUMBER
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE READ_VTK_MESH(INPUT_FILE_NAME,WORLD_REGION,REGION_NUMBER,MESH,GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VARYING_STRING), INTENT(IN) :: INPUT_FILE_NAME
    TYPE(REGION_TYPE), INTENT(IN), POINTER :: WORLD_REGION
    INTEGER(INTG), INTENT(IN) :: REGION_NUMBER
    TYPE(MESH_TYPE), INTENT(INOUT), POINTER :: MESH
    TYPE(FIELD_TYPE), INTENT(INOUT), POINTER :: GEOMETRIC_FIELD
    TYPE(DECOMPOSITION_TYPE), INTENT(INOUT), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    INTEGER(INTG) :: TEXT_LENGTH, NUMBER_OF_CELLS, &
        & NUMBER_OF_ELEMENT_COMPONENTS, NUMBER_OF_PROCESSORS, TEMP_INT, i, j, k, NUMBER_OF_CELL_VALUES, SPACE_DIMENSION=3, &
        & FILE_HANDLE=113
    INTEGER(LINTG), ALLOCATABLE :: CELLS(:), CELL_TYPES(:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    REAL(DP), ALLOCATABLE :: NODES_COORDINATES(:), NODE_COORDINATES(:)
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING, TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3, TEMP_STRING_WORD_4
    character(LEN=MAX_CHARS_PER_LINE) :: TEMP_STRING_CLEAN
    TYPE(VARYING_STRING) :: EXPECTED_LINE
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(BASIS_TYPE), POINTER :: BASIS
    LOGICAL :: EOF, NODES_READ, CELLS_READ, CELL_TYPES_READ
    TYPE(VTK_FIELDS_SET_TYPE) :: FIELD_DATA
    TYPE(VARYING_STRING) :: EVS ! Empty varying string
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("READ_VTK_MESH",ERR,ERROR,*999)
    
    EVS=""
    
    DO i=1,len(TEMP_STRING_CLEAN,1)
      TEMP_STRING_CLEAN(i:i)=" "
    ENDDO
    
    OPEN (FILE_HANDLE,FILE=CHAR(INPUT_FILE_NAME))
    
    ! ----- first line: file format version -----
    EXPECTED_LINE="# vtk DataFile Version 3.0"
    TEMP_STRING=TEMP_STRING_CLEAN
    READ (FILE_HANDLE,'(A)') TEMP_STRING
    IF (TRIM(TEMP_STRING) /= EXPECTED_LINE) THEN
      LOCAL_ERROR='Invalid VTK file version. Only version 3.0 is supported, therefore line 1 '// &
          & 'of your VTK file should be "'//EXPECTED_LINE//'", but found "'//TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    
    ! ----- second line: header -----
    TEMP_STRING=TEMP_STRING_CLEAN
    READ (FILE_HANDLE,'(A)') TEMP_STRING
    
    ! ----- third line: file encoding (ASCII or binary) -----
    TEMP_STRING=TEMP_STRING_CLEAN
    READ (FILE_HANDLE,'(A)') TEMP_STRING
    IF (TRIM(TEMP_STRING) == "ASCII") THEN
      !OK
    ELSEIF (TRIM(TEMP_STRING) == "BINARY") THEN
      LOCAL_ERROR='Invalid file encoding. Only ASCII is supported, but your file has binary type.'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      LOCAL_ERROR='Found invalid keyword for file encoding in line 3, expected "ASCII" or "BINARY", but got "'// &
          & TRIM(TEMP_STRING)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    
    ! ----- fourth line: file type -----
    TEMP_STRING=TEMP_STRING_CLEAN
    READ (FILE_HANDLE,'(A)') TEMP_STRING
    TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
    TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
    READ (TEMP_STRING,*) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2
    IF (TRIM(TEMP_STRING_WORD_1) /= "DATASET") THEN
      LOCAL_ERROR='Expected keyword "DATASET" in line 4, but found "'//TRIM(TEMP_STRING_WORD_1)//'".'
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    SELECT CASE(TRIM(TEMP_STRING_WORD_2))
      CASE("UNSTRUCTURED_GRID")
        TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
        TEMP_STRING_WORD_2=TEMP_STRING_CLEAN
        TEMP_STRING_WORD_3=TEMP_STRING_CLEAN
        READ (TEMP_STRING,*,end=881) TEMP_STRING_WORD_1, TEMP_STRING_WORD_2, TEMP_STRING_WORD_3
        !Still here (no jump to 881) => TEMP_STRING_WORD_3 was found => Error
        LOCAL_ERROR="UNSTRUCTURED_GRID type doesn't take any additional parameters, "// &
            & 'but found "'//TRIM(TEMP_STRING_WORD_3)//'".'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
881     CONTINUE !everything OK, TEMP_STRING_WORD_3 was not found
        EOF=.False.
        NODES_READ=.False.
        CELLS_READ=.False.
        CELL_TYPES_READ=.False.
        ALLOCATE(FIELD_DATA%POINT_ASSOCIATED_FIELDS(0))
        ALLOCATE(FIELD_DATA%CELL_ASSOCIATED_FIELDS(0))
        DO WHILE(.not. EOF)
          TEMP_STRING=TEMP_STRING_CLEAN
          READ (FILE_HANDLE,'(A)',END=883) TEMP_STRING
          TEMP_STRING_WORD_1=TEMP_STRING_CLEAN
          READ (TEMP_STRING,*,end=882) TEMP_STRING_WORD_1
          SELECT CASE(TRIM(TEMP_STRING_WORD_1))
            CASE("POINTS")
              CALL READ_VTK_POINTS(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,NODES_COORDINATES,ERR,ERROR,*999)
              NODES_READ=.True.
            CASE("CELLS")
              CALL READ_VTK_CELLS(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,CELLS,ERR,ERROR,*999)
              CELLS_READ=.True.
            CASE("CELL_TYPES")
              CALL READ_VTK_CELL_TYPES(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,CELL_TYPES,ERR,ERROR,*999)
              CELL_TYPES_READ=.True.
            CASE("POINT_DATA")
              CALL READ_VTK_CELL_AND_POINT_DATA(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,FIELD_DATA,ERR,ERROR,*999)
              EOF=.True.
            CASE("CELL_DATA")
              CALL READ_VTK_CELL_AND_POINT_DATA(FILE_HANDLE,TEMP_STRING,TEMP_STRING_CLEAN,FIELD_DATA,ERR,ERROR,*999)
              EOF=.True.
            CASE DEFAULT
              LOCAL_ERROR='Found invalid keyword "'//TRIM(TEMP_STRING_WORD_1)// &
                  & '" for section name. Expected POINTS, CELLS, CELL_TYPES, POINT_DATA or CELL_DATA'
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
882       CONTINUE ! empty line
        ENDDO
883     CONTINUE ! EOF
      CASE("STRUCTURED_POINTS", "STRUCTURED_GRID", "POLYDATA", "RECTILINEAR_GRID", "FIELD")
        LOCAL_ERROR='Dataset type '//TRIM(TEMP_STRING_WORD_2)//' not supported.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR='Found invalid keyword "'//TRIM(TEMP_STRING)//'" for dataset type in line 4.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    
    ! Create coordinate system
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_CREATE_START(REGION_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)
    
    ! Create region and assign coordinate system to it
    NULLIFY(REGION)
    CALL REGION_CREATE_START(REGION_NUMBER,WORLD_REGION,REGION,ERR,ERROR,*999)
    CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)
    
    ! Create linear basis
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START(REGION_NUMBER,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_SET(BASIS,BASIS_SIMPLEX_TYPE,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
        & BASIS_LINEAR_SIMPLEX_INTERPOLATION/),ERR,ERROR,*999)
    CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
    
    !Create mesh
    NULLIFY(MESH)
    CALL MESH_CREATE_START(REGION_NUMBER,REGION,SPACE_DIMENSION,MESH,ERR,ERROR,*999)
    
    ! Create nodes
    NULLIFY(NODES)
    IF (.NOT. NODES_READ) THEN
        LOCAL_ERROR='No nodes found in VTK file.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    CALL NODES_CREATE_START(REGION,SIZE(NODES_COORDINATES,1)/SPACE_DIMENSION,NODES,ERR,ERROR,*999)
    CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
    
    ! Create elements
    IF (.NOT. CELLS_READ) THEN
        LOCAL_ERROR='No cells found in VTK file.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    IF (.NOT. CELL_TYPES_READ) THEN
        LOCAL_ERROR='No cell types found in VTK file.'
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,1,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,size(CELL_TYPES,1),ERR,ERROR,*999)
    NULLIFY(ELEMENTS)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,ELEMENTS,ERR,ERROR,*999)
    j=1
    DO i=1,SIZE(CELL_TYPES,1)
      SELECT CASE(CELL_TYPES(i))
        CASE(VTK_EMPTY_CELL)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_EMPTY_CELL",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_VERTEX)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_VERTEX",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_POLY_VERTEX)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_POLY_VERTEX",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_LINE)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_LINE",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_POLY_LINE)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_POLY_LINE",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_TRIANGLE)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_TRIANGLE",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_TRIANGLE_STRIP)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_TRIANGLE_STRIP",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_POLYGON)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_POLYGON",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_PIXEL)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_PIXEL",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_QUAD)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_QUAD",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_TETRA)
          CALL READ_VTK_MESH_CHECK_CELL_NODE_NUMBER(EVS//"VTK_TETRA",4,CELLS(j),i,ERR,ERROR,*999)
          ALLOCATE(ELEMENT_NODES(CELLS(j)),STAT=ERR)
          DO k=j+1,j+CELLS(j)
            ELEMENT_NODES(k-j)=CELLS(k)+1
          ENDDO
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(i,ELEMENTS,ELEMENT_NODES,Err,ERROR,*999)
          DEALLOCATE(ELEMENT_NODES)
          j=j+CELLS(j)+1
        CASE(VTK_VOXEL)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_VOXEL",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_HEXAHEDRON)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_HEXAHEDRON",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_WEDGE)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_WEDGE",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_PYRAMID)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_PYRAMID",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_PENTAGONAL_PRISM)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_PENTAGONAL_PRISM",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE(VTK_HEXAGONAL_PRISM)
          CALL READ_VTK_MESH_ERROR_CELL_TYPE(EVS//"VTK_HEXAGONAL_PRISM",CELL_TYPES(i),i,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR='Cell type '//TRIM(NUMBER_TO_VSTRING(CELL_TYPES(i),"*",ERR,ERROR))// &
              & ' invalid (found in cell '//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//')!'
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ENDDO
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*999)
    
    CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999)
    
    !Calculate decomposition
    NULLIFY(DECOMPOSITION)
    NUMBER_OF_PROCESSORS = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    CALL DECOMPOSITION_CREATE_START(REGION_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_PROCESSORS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*999)

    !Create a field to put the geometry
    NULLIFY(GEOMETRIC_FIELD)
    CALL FIELD_CREATE_START(REGION_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    CALL FIELD_TYPE_SET(GEOMETRIC_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_VARIABLES_SET(GEOMETRIC_FIELD,1,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_COMPONENTS_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,SPACE_DIMENSION,Err,ERROR,*999)
    DO i=1,SPACE_DIMENSION
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,i,1,ERR,ERROR,*999)
    ENDDO
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)

    !Set node positions
    ALLOCATE(NODE_COORDINATES(SPACE_DIMENSION),STAT=ERR)
    DO i=1,size(NODES_COORDINATES,1)/SPACE_DIMENSION
      NODE_COORDINATES=NODES_COORDINATES(SPACE_DIMENSION*(i-1)+1:SPACE_DIMENSION*i)
      DO j=1,SPACE_DIMENSION
        CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,i,j, &
          & NODE_COORDINATES(j),ERR,ERROR,*999)
      ENDDO
    ENDDO
    DEALLOCATE(NODE_COORDINATES)
    
    CLOSE(FILE_HANDLE)
    
    CALL EXITS("READ_VTK_MESH")
    RETURN
999 CALL ERRORS("READ_VTK_MESH",ERR,ERROR)
    CALL EXITS("READ_VTK_MESH")
    RETURN 1
    
  END SUBROUTINE READ_VTK_MESH
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE WRITE_VTK_MESH(OUTPUT_FILE_NAME,MESH,GEOMETRIC_FIELD,OUTPUT_FIELDS,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VARYING_STRING), INTENT(IN) :: OUTPUT_FILE_NAME
    TYPE(MESH_TYPE), INTENT(IN), POINTER :: MESH
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: GEOMETRIC_FIELD
    TYPE(FIELD_PTR_TYPE), INTENT(IN), ALLOCATABLE :: OUTPUT_FIELDS(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS, NUMBER_OF_NODES, node_idx, dim_idx, local_ny, ne, &
            & NUMBER_OF_NODES_PER_ELEMENT, i
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: GEOMETRIC_VARIABLE,OUTPUT_FIELD_VARIABLE
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:), OUTPUT_FIELD_PARAMETERS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS
    REAL(DP) :: value
    
    CALL ENTERS("WRITE_VTK_MESH",ERR,ERROR,*999)
    
    OPEN (12,FILE=CHAR(OUTPUT_FILE_NAME // ".vtk"))
    
    WRITE(12,'(A)') "# vtk DataFile Version 3.0"
    WRITE(12,'(A)') "vtk output"
    WRITE(12,'(A)') "ASCII"
    WRITE(12,'(A)') "DATASET UNSTRUCTURED_GRID"
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    NULLIFY(GEOMETRIC_VARIABLE)
    CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
            & ERR,ERROR,*999)
    
    NUMBER_OF_NODES=GEOMETRIC_VARIABLE%COMPONENTS(1)%DOMAIN%TOPOLOGY%NODES%NUMBER_OF_NODES
    WRITE(12,'(A,I8,A6)') "POINTS",NUMBER_OF_NODES,"float"
    
    DO node_idx=1,NUMBER_OF_NODES
      DO dim_idx=1,NUMBER_OF_DIMENSIONS
        local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
            & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
        WRITE(12,*) GEOMETRIC_PARAMETERS(local_ny)
      ENDDO
    ENDDO
    
    ELEMENTS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS
    NUMBER_OF_NODES_PER_ELEMENT=size(ELEMENTS%ELEMENTS(1)%ELEMENT_NODES,1)
    WRITE(12,'(A,I8,I8)') "CELLS ",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,ELEMENTS% &
            & TOTAL_NUMBER_OF_ELEMENTS*(NUMBER_OF_NODES_PER_ELEMENT+1)
    DO ne=1,ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,*) NUMBER_OF_NODES_PER_ELEMENT, &
            & ((ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(i)-1),i=1,size(ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES,1))
    ENDDO
    
    WRITE(12,'(A,I8)') "CELL_TYPES",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
    DO ne=1,ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,'(A)') "10"
    ENDDO
    
    WRITE(12,'(A,I8)') "CELL_DATA",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
    WRITE(12,'(A,I8)') "POINT_DATA",NUMBER_OF_NODES
    
    DO i=1,SIZE(OUTPUT_FIELDS, 1)
      CALL FIELD_NUMBER_OF_COMPONENTS_GET(OUTPUT_FIELDS(i)%PTR,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
      !NULLIFY(OUTPUT_FIELD_VARIABLE)
      !CALL FIELD_VARIABLE_GET(OUTPUT_FIELDS(i)%PTR,FIELD_U_VARIABLE_TYPE,OUTPUT_FIELD_VARIABLE,ERR,ERROR,*999)
      !CALL FIELD_PARAMETER_SET_DATA_GET(OUTPUT_FIELDS(i)%PTR,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      !        & OUTPUT_FIELD_PARAMETERS,ERR,ERROR,*999)
      
      WRITE(12,'(A,A)') "FIELD number"," 1"
      WRITE(12,'(I3,I3,I8,A6)') i,1,NUMBER_OF_NODES,"float"
      DO node_idx=1,NUMBER_OF_NODES
        DO dim_idx=1,NUMBER_OF_DIMENSIONS
          CALL FIELD_PARAMETER_SET_GET_NODE(OUTPUT_FIELDS(i)%PTR,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,1,1,node_idx,1,value,ERR,ERROR,*999)
          !local_ny=OUTPUT_FIELD_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
          !    & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
          WRITE(12,*) value
          !WRITE(12,*) OUTPUT_FIELD_PARAMETERS(local_ny)
        ENDDO
      ENDDO
    ENDDO
    
    !      export FIELD information
    !WRITE(12,'(A,A)')"FIELD number"," 1"
    !WRITE(12,'(A,I3,I8,A6)')OUTPUT_FILE_FIELD_TITLE,1,TOTAL_NUMBER_OF_NODES,"float"
    !DO I=1,TOTAL_NUMBER_OF_NODES
    !WRITE(12,'(F15.10)') SEED_VALUE(I)
    !ENDDO
    
    !      export VECTORS information
    !WRITE(12,'(A,A,A6)') "VECTORS ","fiber_vector","float"
    !DO I=1,TOTAL_NUMBER_OF_NODES
    !WRITE(12,'(3F8.5)') (CONDUCTIVITY_TENSOR(I,J),J=1,3)
    !ENDDO
    
    CLOSE(12)
    
    CALL EXITS("WRITE_VTK_MESH")
    RETURN
999 CALL ERRORS("WRITE_VTK_MESH",ERR,ERROR)
    CALL EXITS("WRITE_VTK_MESH")
    RETURN 1
    
  END SUBROUTINE WRITE_VTK_MESH
  
  
  
END MODULE MESH_IO_ROUTINES


/*
 * \file
 */

#ifndef __PDM_PREDICATE_H__
#define __PDM_PREDICATE_H__

/*****************************************************************************/
/*                                                                           */
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
/*  and Fast Robust Geometric Predicates                                     */
/*  (predicates.c)                                                           */
/*                                                                           */
/*  May 18, 1996                                                             */
/*                                                                           */
/*  Placed in the public domain by                                           */
/*  Jonathan Richard Shewchuk                                                */
/*  School of Computer Science                                               */
/*  Carnegie Mellon University                                               */
/*  5000 Forbes Avenue                                                       */
/*  Pittsburgh, Pennsylvania  15213-3891                                     */
/*  jrs@cs.cmu.edu                                                           */
/*                                                                           */
/*  This file contains C implementation of algorithms for exact addition     */
/*    and multiplication of floating-point numbers, and predicates for       */
/*    robustly performing the orientation and incircle tests used in         */
/*    computational geometry.  The algorithms and underlying theory are      */
/*    described in Jonathan Richard Shewchuk.  "Adaptive Precision Floating- */
/*    Point Arithmetic and Fast Robust Geometric Predicates."  Technical     */
/*    Report CMU-CS-96-140, School of Computer Science, Carnegie Mellon      */
/*    University, Pittsburgh, Pennsylvania, May 1996.  (Submitted to         */
/*    Discrete & Computational Geometry.)                                    */
/*                                                                           */
/*  This file, the paper listed above, and other information are available   */
/*    from the Web page http://www.cs.cmu.edu/~quake/robust.html .           */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Using this code:                                                         */
/*                                                                           */
/*  First, read the short or long version of the paper (from the Web page    */
/*    above).                                                                */
/*                                                                           */
/*  Be sure to call exactinit() once, before calling any of the arithmetic   */
/*    functions or geometric predicates.  Also be sure to turn on the        */
/*    optimizer when compiling this file.                                    */
/*                                                                           */
/*                                                                           */
/*  Several geometric predicates are defined.  Their parameters are all      */
/*    points.  Each point is an array of two or three floating-point         */
/*    numbers.  The geometric predicates, described in the papers, are       */
/*                                                                           */
/*    orient2d(pa, pb, pc)                                                   */
/*    orient2dfast(pa, pb, pc)                                               */
/*    orient3d(pa, pb, pc, pd)                                               */
/*    orient3dfast(pa, pb, pc, pd)                                           */
/*    incircle(pa, pb, pc, pd)                                               */
/*    incirclefast(pa, pb, pc, pd)                                           */
/*    insphere(pa, pb, pc, pd, pe)                                           */
/*    inspherefast(pa, pb, pc, pd, pe)                                       */
/*                                                                           */
/*  Those with suffix "fast" are approximate, non-robust versions.  Those    */
/*    without the suffix are adaptive precision, robust versions.  There     */
/*    are also versions with the suffices "exact" and "slow", which are      */
/*    non-adaptive, exact arithmetic versions, which I use only for timings  */
/*    in my arithmetic papers.                                               */
/*                                                                           */
/*                                                                           */
/*  An expansion is represented by an array of floating-point numbers,       */
/*    sorted from smallest to largest magnitude (possibly with interspersed  */
/*    zeros).  The length of each expansion is stored as a separate integer, */
/*    and each arithmetic function returns an integer which is the length    */
/*    of the expansion it created.                                           */
/*                                                                           */
/*  Several arithmetic functions are defined.  Their parameters are          */
/*                                                                           */
/*    e, f           Input expansions                                        */
/*    elen, flen     Lengths of input expansions (must be >= 1)              */
/*    h              Output expansion                                        */
/*    b              Input scalar                                            */
/*                                                                           */
/*  The arithmetic functions are                                             */
/*                                                                           */
/*    grow_expansion(elen, e, b, h)                                          */
/*    grow_expansion_zeroelim(elen, e, b, h)                                 */
/*    expansion_sum(elen, e, flen, f, h)                                     */
/*    expansion_sum_zeroelim1(elen, e, flen, f, h)                           */
/*    expansion_sum_zeroelim2(elen, e, flen, f, h)                           */
/*    fast_expansion_sum(elen, e, flen, f, h)                                */
/*    fast_expansion_sum_zeroelim(elen, e, flen, f, h)                       */
/*    linear_expansion_sum(elen, e, flen, f, h)                              */
/*    linear_expansion_sum_zeroelim(elen, e, flen, f, h)                     */
/*    scale_expansion(elen, e, b, h)                                         */
/*    scale_expansion_zeroelim(elen, e, b, h)                                */
/*    compress(elen, e, h)                                                   */
/*                                                                           */
/*  All of these are described in the long version of the paper; some are    */
/*    described in the short version.  All return an integer that is the     */
/*    length of h.  Those with suffix _zeroelim perform zero elimination,    */
/*    and are recommended over their counterparts.  The procedure            */
/*    fast_expansion_sum_zeroelim() (or linear_expansion_sum_zeroelim() on   */
/*    processors that do not use the round-to-even tiebreaking rule) is      */
/*    recommended over expansion_sum_zeroelim().  Each procedure has a       */
/*    little note next to it (in the code below) that tells you whether or   */
/*    not the output expansion may be the same array as one of the input     */
/*    expansions.                                                            */
/*                                                                           */
/*                                                                           */
/*  If you look around below, you'll also find macros for a bunch of         */
/*    simple unrolled arithmetic operations, and procedures for printing     */
/*    expansions (commented out because they don't work with all C           */
/*    compilers) and for generating random floating-point numbers whose      */
/*    significand bits are all random.  Most of the macros have undocumented */
/*    requirements that certain of their parameters should not be the same   */
/*    variable; for safety, better to make sure all the parameters are       */
/*    distinct variables.  Feel free to send email to jrs@cs.cmu.edu if you  */
/*    have questions.                                                        */
/*                                                                           */
/*****************************************************************************/



/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 *  from Triangle (https://github.com/libigl/triangle)
 */


#define REAL double
#define vertex double *

#define INEXACT /* Nothing */
/* #define INEXACT volatile */
void PDM_predicate_exactinit (void);


/**
 *
 * \brief Determine if a point D is inside/out/on the circumcircle of triangle ABC
 *
 * This function assumes that all points lie in the z=0 plane.
 *
 * \param [in]   pa   Point A
 * \param [in]   pb   Point B
 * \param [in]   pc   Point C
 * \param [in]   pd   Point D
 *
 * \return Strictly positive value if D is inside the circumdisk
 *         Strictly negative value if D is outside the circumdisk
 *         Zero if D is on the circumcircle
 *
 */

REAL PDM_predicate_incircle
(
 vertex pa,
 vertex pb,
 vertex pc,
 vertex pd
 );


/**
 *
 * \brief Determine if a point E is inside/out/on the circumsphere of tetrahedron ABCD
 *
 * \param [in]   pa   Point A
 * \param [in]   pb   Point B
 * \param [in]   pc   Point C
 * \param [in]   pd   Point D
 * \param [in]   pe   Point E
 *
 * \return Strictly positive value if E is inside the circumball
 *         Strictly negative value if E is outside the circumball
 *         Zero if E is on the circumsphere
 *
 */

REAL PDM_predicate_insphere
(
 REAL *pa,
 REAL *pb,
 REAL *pc,
 REAL *pd,
 REAL *pe
 );


/**
 *
 * \brief Determine the orientation of triangle ABC.
 *
 * This function assumes that all points lie in the z=0 plane.
 *
 * \param [in]   pa   Point A
 * \param [in]   pb   Point B
 * \param [in]   pc   Point C
 *
 * \return Twice the signed area of triangle ABC
 *
 */

REAL PDM_predicate_orient2d
(
 REAL *pa,
 REAL *pb,
 REAL *pc
 );



/**
 *
 * \brief Determine the orientation of tetrahedron ABCD.
 * *
 * \param [in]   pa   Point A
 * \param [in]   pb   Point B
 * \param [in]   pc   Point C
 * \param [in]   pd   Point D
 *
 * \return Twice the signed volume of tetrahedron ABCD
 *
 */

REAL PDM_predicate_orient3d
(
 REAL *pa,
 REAL *pb,
 REAL *pc,
 REAL *pd
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PREDICATE_H__ */

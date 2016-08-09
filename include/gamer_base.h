/**
 * @defgroup global_gamer global_gamer class
 * @brief    Global group for gamer
 */

/**
 *  @file       gamer_base.h
 *  @ingroup    global_gamer
 *  @brief      The base (or foundation) header for GAMER.
 *  @authors    Michael Holst
 *  @note       This header sets things up correctly for using ISO/ANSI-C.
 *              The following macros affect the behavior of the header:     
    @verbatim
     Inlining for speed:  (Normal C functions if VINLINE_XXX not defined.)    
     ------------------
     -DVINLINE_GAMER : Enables macro time-critical funcs in GAMER.
    @endverbatim                                                             
 *
 *  @version    $Id: gamer_base.h,v 1.11 2010/08/12 05:43:10 fetk Exp $
 *  
 *  @attention
 *  @verbatim
 *
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 1994-- Michael Holst and Zeyun Yu
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
 */

#ifndef _GAMER_BASE_H_
#define _GAMER_BASE_H_

#include <maloc/maloc.h>

/** @brief Triangle-specific macros */
#define TRILIBRARY
/** @brief Triangle-specific macros */
#define REAL double
/** @brief Triangle-specific macros */
#define VOID void
/** @brief Triangle-specific macros */
#define ANSI_DECLARATORS


/** @brief Tetgen-specific macros */
#define TETLIBRARY
/* #define REAL double */  /* This is common with Triangle */

/*
 * ***************************************************************************
 * Biom-specific macros
 * ***************************************************************************
 */

#endif /* _GAMER_BASE_H_ */


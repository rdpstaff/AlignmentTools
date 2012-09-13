/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.msu.cme.rdp.alignment;

/**
 *
 * @author fishjord
 */
public enum AlignmentMode {
    /**Local alignment only*/local,
    /**Align seqi and seqj from seqi's start and end points*/glocal,
    /**Align sequence with no leading/trailing gap penalty*/overlap,
    /**Align sequences with no leading/training gap penalty, only return the overlapping region*/overlap_trim,
    /**Global alignment*/global
}

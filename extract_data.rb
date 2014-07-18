#!/usr/bin/env ruby

# Extract 2D slices or 3D subsets from binary output of ETD Stepper.
# Author: Matt Pulver <matt@blue-science.org>
# Date: 2014 July
#
# Usage: zcat frames/frame000.000000.gz | ./extract_data.rb space > frame00.bin
# where space = xy, yz, xz, or xyz:
#   xy  : z=0 plane
#   yz  : x=0 plane
#   xz  : y=0 plane
#   xyz : subset of points determined by filter condition below
# It is assumed that each spatial dimension is centered at 0.

# STDIN: Output from ETD Stepper
# STDOUT: Gnuplot matrix of phi_Y values in 2D or 3D.

T, RANK = STDIN.binmode.read(12).unpack 'df'
abort 'This script is only for 3-dimensional input.' if RANK != 3
L, DIMS = STDIN.read(RANK*8).unpack('ff'*RANK).each_slice(2).to_a.transpose
N = DIMS.inject 1, &:*
STDERR.puts "N = #{N}, T = #{T}"

space = ARGV[0] || 'xy'
abort "Usage: zcat frames/frameMMM.MMMMMM.gz | #{__FILE__} [xy|yz|xz|xyz]" if !/^(xy|yz|xz|xyz)$/.match(space)

STDOUT.binmode

# Output data for the plane x = 0
def space_yz
    # Print first row. These are the z-values.
    STDOUT.print ([DIMS[2]]+(0...DIMS[2]).map { |iz| (iz-DIMS[2]/2)*L[2]/DIMS[2] }).pack('f*')

    # The plane x = 0 corresponds with ix = DIMS[0]/2 (N/DIMS[0]*(DIMS[0]/2) ... 
    start = N/DIMS[0]*(DIMS[0].to_i/2) # integer division
    STDIN.read start*12

    (0...N/DIMS[0]).each do |n|
        iy, iz = n.divmod DIMS[2].to_i
        pg, px, py = STDIN.read(12).unpack('fff')
        STDOUT.print [(iy-DIMS[1]/2)*L[1]/DIMS[1]].pack('f') if iz.zero?
        STDOUT.print [py].pack('f')
    end
end

# Output data for the plane z = 0
def space_xy
    # Print first row. These are the y-values.
    STDOUT.print ([DIMS[1]]+(0...DIMS[1]).map { |iy| (iy-DIMS[1]/2)*L[1]/DIMS[1] }).pack('f*')

    start = 12*(DIMS[2].to_i/2-1) # integer division
    STDIN.read start
    skip = 12*(DIMS[2]-1)

    (0...N/DIMS[2]).each do |n|
        ix, iy = n.divmod DIMS[1].to_i
        STDIN.read skip if !n.zero? # skip to next value for which z=0
        pg, px, py = STDIN.read(12).unpack('fff')
        STDOUT.print [(ix-DIMS[0]/2)*L[0]/DIMS[0]].pack('f') if iy.zero?
        STDOUT.print [py].pack('f')
    end
end

# Output data for the plane y = 0
def space_xz
    # Print first row. These are the z-values.
    STDOUT.print ([DIMS[2]]+(0...DIMS[2]).map { |iz| (iz-DIMS[2]/2)*L[2]/DIMS[2] }).pack('f*')

    start = 12*(DIMS[1].to_i/2)*DIMS[2] # integer division
    STDIN.read start
    skip = 12*(DIMS[1]-1)*DIMS[2]

    (0...N/DIMS[1]).each do |n|
        ix, iz = n.divmod DIMS[2].to_i
        STDIN.read skip if !ix.zero? && iz.zero? # skip to next value for which y=0
        pg, px, py = STDIN.read(12).unpack('fff')
        STDOUT.print [(ix-DIMS[0]/2)*L[0]/DIMS[0]].pack('f') if iz.zero?
        STDOUT.print [py].pack('f')
    end
end

# Output 3D data that matches the below filter condition
def space_xyz
    (0...N).each do |n|
        pg, px, py = STDIN.read(12).unpack('fff')
        next if py < 0.5 # This is the filter condition
        iy, iz =  n.divmod DIMS[2].to_i
        ix, iy = iy.divmod DIMS[1].to_i
        x, y, z = [ix,iy,iz].map.with_index { |i,j| (i-DIMS[j]/2)*L[j]/DIMS[j] }
        STDOUT.print [x,y,z,py].pack('f*')
    end
end

send "space_#{space}"

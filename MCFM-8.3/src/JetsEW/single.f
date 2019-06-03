C ---------------------------------------------------------------------- C
C     Born for quark induced dijets processes. It differentiates         C
C     from channel to channel, namely, each piece of the function        C
C     returns a matrix element of one channel interferes with another,   C
C     which is similar with the functions in "small.f".                  C
C ---------------------------------------------------------------------- C

      function sxs(ss,tt,uu)
C     function sxs denotes s-channel interferes with s-channel, and similar
C     notations are used in the following little functions
C     the s-channel process is defined as qa_ii_jj
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: ss,tt,uu,ssq,tsq,usq,sxs

      ssq = ss**2
      tsq = tt**2
      usq = uu**2

!      sxs = 2._dp*V*(tsq + usq)/ssq
C --- strip off the color structure
      sxs = + V*(tsq + usq)/ssq

      end function sxs



      function sxt(ss,tt,uu)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: ss,tt,uu,usq,sxt

      usq = uu**2

!      sxt = -2._dp*V/xn*usq/(ss*tt)
C --- strip off the color structure; served as mixed born instead of qcd born
C --- refer to routine 'dijet_qqb_born_mix' in 'dijet_qqb_mix.f'
c      sxt = V*usq/(ss*tt)
c      sxt =  - V*uu**2/ss/tt
C --- add '-' due to single fermion-loop in s- and t-channel interference
      sxt =  + V*uu**2/ss/tt

      end function sxt


C -- list the crossing symmetries at tree level
C     1. qa_ii_jj = sxs(ss,tt,uu) ~ (tsq + usq)/ssq
C     2. aq_ii_jj = sxs(ss,uu,tt) ~ (tsq + usq)/ssq
C     3. qq_ij_ij = sxs(tt,uu,ss) ~ (ssq + usq)/tsq
C     4. aa_ij_ij = sxs(tt,uu,ss) ~ (ssq + usq)/tsq
C     5. qa_ij_ij = sxs(tt,ss,uu) ~ (ssq + usq)/tsq
C     6. aq_ij_ij = sxs(tt,ss,uu) ~ (ssq + usq)/tsq
C     7. qa_ii_ii = sxs(ss,tt,uu) + sxs(tt,ss,uu) + 2*sxt(ss,tt,uu)
C                 ~ smallb(uu,ss,tt)
C     8. aq_ii_ii = sxs(ss,tt,uu) + sxs(tt,uu,ss) + 2*sxt(ss,tt,uu)
C                 ~ smallb(uu,ss,tt)
C     9. qq_ii_ii = sxs(tt,uu,ss) + sxs(uu,tt,ss) + 2*sxt(tt,uu,ss)
C                 ~ smallb(ss,tt,uu)
C     10.aa_ii_ii = sxs(tt,uu,ss) + sxs(uu,tt,ss) + 2*sxt(tt,uu,ss)
C                 ~ smallb(ss,tt,uu)

C     since we know the vertex correction is proportional to the Born, 
C     and the correction is contained in the so-called form factor, we
C     can use the Born and the form factor information to get the vertex
C     correction.

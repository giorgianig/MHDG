syms xx yy xm ym D mu csii csie kpari kpare a cx cy sx sy cx2 cy2 sx2 sy2 r r2
syms epn Mref pcourr tie cc ss

fortran( ((10*a*cx2*cy*(ym-yy)*(cy-5*sy+cy*ss))/(3*((xm^2-2*ym*yy-...
            2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^(1/2))-(5*a*cy*sx*(cy*sx-10)*(ss+2)*...
            (ym-yy))/(3*((xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^(1/2))+...
            (5*cc*(cy*sx-10)*(ss+2)*(ym-yy)*(xm^2-2*ym*yy-xm*xx+ym^2+yy^2))/...
            (3*xx^3*((xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^(3/2)))/xx+...
            (4*3^(1-epn)*Mref*a*kpare*ym^4*cc*(-(2*cy*sx-20)/Mref)^epn+4*3^(1-epn)*...
            Mref*a*kpare*yy^4*cc*(-(2*cy*sx-20)/Mref)^epn-2*3^(1-epn)*Mref*...
            a^2*kpare*xx^3*cy*sx*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*...
            xx+2*xx^2+ym^2+yy^2)-16*3^(1-epn)*Mref*a*kpare*ym*yy^3*cc*(-(2*cy*sx-20)/...
            Mref)^epn-16*3^(1-epn)*Mref*a*kpare*ym^3*yy*cc*(-(2*cy*sx-20)/...
            Mref)^epn+4*3^(1-epn)*Mref*a*kpare*xm*ym^3*ss*(-(2*cy*sx-20)/Mref)^epn+...
            4*3^(1-epn)*Mref*a*kpare*xm^3*ym*ss*(-(2*cy*sx-20)/Mref)^epn-4*...
            3^(1-epn)*Mref*a*kpare*xx*ym^3*ss*(-(2*cy*sx-20)/Mref)^epn-4*3^...
            (1-epn)*Mref*a*kpare*xx^3*ym*ss*(-(2*cy*sx-20)/Mref)^epn-4*3^...
            (1-epn)*Mref*a*kpare*xm*yy^3*ss*(-(2*cy*sx-20)/Mref)^epn-4*3^...
            (1-epn)*Mref*a*kpare*xm^3*yy*ss*(-(2*cy*sx-20)/Mref)^epn+4*3^...
            (1-epn)*Mref*a*kpare*xx*yy^3*ss*(-(2*cy*sx-20)/Mref)^epn+4*3^...
            (1-epn)*Mref*a*kpare*xx^3*yy*ss*(-(2*cy*sx-20)/Mref)^epn-2*3^...
            (1-epn)*Mref*a*kpare*xx^2*cc*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*...
            ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)-2*3^(1-epn)*Mref*a*kpare*ym^2*...
            cc*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+...
            yy^2)-2*3^(1-epn)*Mref*a*kpare*yy^2*cc*(-(2*cy*sx-20)/Mref)^epn*...
            (xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)+4*3^(1-epn)*Mref*a*kpare*...
            xm^2*ym^2*cc*(-(2*cy*sx-20)/Mref)^epn+4*3^(1-epn)*Mref*a*kpare*...
            xx^2*ym^2*cc*(-(2*cy*sx-20)/Mref)^epn+4*3^(1-epn)*Mref*a*kpare*...
            xm^2*yy^2*cc*(-(2*cy*sx-20)/Mref)^epn+4*3^(1-epn)*Mref*a*kpare*...
            xx^2*yy^2*cc*(-(2*cy*sx-20)/Mref)^epn+8*3^(2-epn)*Mref*a*kpare*...
            ym^2*yy^2*cc*(-(2*cy*sx-20)/Mref)^epn+4*3^(1-epn)*Mref*a^2*kpare*...
            xm*xx^2*cy*sx*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*...
            xx^2+ym^2+yy^2)-2*3^(1-epn)*Mref*a^2*kpare*xm^2*xx*cy*sx*(-(2*cy*...
            sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)-2*3^...
            (1-epn)*Mref*a^2*kpare*xx*ym^2*cy*sx*(-(2*cy*sx-20)/Mref)^epn*...
            (xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)-4*3^(1-epn)*Mref*a^2*...
            kpare*xx^2*ym*cx*sy*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*...
            xx+2*xx^2+ym^2+yy^2)-2*3^(1-epn)*Mref*a^2*kpare*xx*yy^2*cy*sx*...
            (-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)+...
            4*3^(1-epn)*Mref*a^2*kpare*xx^2*yy*cx*sy*(-(2*cy*sx-20)/Mref)^...
            epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)+2*3^(1-epn)*Mref*a*...
            kpare*xm*xx*cc*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*...
            xx^2+ym^2+yy^2)+4*3^(1-epn)*Mref*a*kpare*ym*yy*cc*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)-8*3^(1-epn)*...
            Mref*a*kpare*xm*xx*ym^2*cc*(-(2*cy*sx-20)/Mref)^epn-8*3^(1-epn)*...
            Mref*a*kpare*xm*xx*yy^2*cc*(-(2*cy*sx-20)/Mref)^epn-8*3^(1-epn)*...
            Mref*a*kpare*xm^2*ym*yy*cc*(-(2*cy*sx-20)/Mref)^epn-8*3^(1-epn)*...
            Mref*a*kpare*xx^2*ym*yy*cc*(-(2*cy*sx-20)/Mref)^epn-2*3^(1-epn)*...
            Mref*a*kpare*xm*ym*ss*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*...
            xm*xx+2*xx^2+ym^2+yy^2)+2*3^(1-epn)*Mref*a*kpare*xm*yy*ss*(-(2*...
            cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)+4*...
            3^(2-epn)*Mref*a*kpare*xm*xx^2*ym*ss*(-(2*cy*sx-20)/Mref)^epn-4*...
            3^(2-epn)*Mref*a*kpare*xm^2*xx*ym*ss*(-(2*cy*sx-20)/Mref)^epn-4*...
            3^(2-epn)*Mref*a*kpare*xm*xx^2*yy*ss*(-(2*cy*sx-20)/Mref)^epn+4*...
            3^(2-epn)*Mref*a*kpare*xm^2*xx*yy*ss*(-(2*cy*sx-20)/Mref)^epn+4*...
            3^(2-epn)*Mref*a*kpare*xm*ym*yy^2*ss*(-(2*cy*sx-20)/Mref)^epn-4*...
            3^(2-epn)*Mref*a*kpare*xm*ym^2*yy*ss*(-(2*cy*sx-20)/Mref)^epn-4*...
            3^(2-epn)*Mref*a*kpare*xx*ym*yy^2*ss*(-(2*cy*sx-20)/Mref)^epn+4*...
            3^(2-epn)*Mref*a*kpare*xx*ym^2*yy*ss*(-(2*cy*sx-20)/Mref)^epn+4*...
            3^(1-epn)*Mref*a^2*kpare*xm*xx*ym*cx*sy*(-(2*cy*sx-20)/Mref)^epn*...
            (xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)-4*3^(1-epn)*Mref*a^2*kpare*...
            xm*xx*yy*cx*sy*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*...
            xx^2+ym^2+yy^2)+4*3^(1-epn)*Mref*a^2*kpare*xx*ym*yy*cy*sx*(-(2*cy*...
            sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)+(4*3^(1-epn)*...
            Mref*a^2*epn*kpare*xx^3*sx2*sy2*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*...
            ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)+16*3^(1-epn)*Mref*a*...
            kpare*xm*xx*ym*yy*cc*(-(2*cy*sx-20)/Mref)^epn+(4*3^(1-epn)*Mref*a^2*...
            epn*kpare*xx*ym^2*cx2*cy2*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*...
            xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)+(4*3^(1-epn)*Mref*a^2*epn*kpare*xx*...
            yy^2*cx2*cy2*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+...
            ym^2+yy^2))/(2*cy*sx-20)-(8*3^(1-epn)*Mref*a^2*epn*kpare*xm*xx^2*sx2*...
            sy2*(-(2*cy*sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/...
            (2*cy*sx-20)+(4*3^(1-epn)*Mref*a^2*epn*kpare*xm^2*xx*sx2*sy2*(-(2*cy*...
            sx-20)/Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)-...
            (8*3^(1-epn)*Mref*a^2*epn*kpare*xx*ym*yy*cx2*cy2*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)-...
            (8*3^(1-epn)*Mref*a^2*epn*kpare*xx^2*ym*cc*ss*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)+...
            (8*3^(1-epn)*Mref*a^2*epn*kpare*xx^2*yy*cc*ss*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)+...
            (8*3^(1-epn)*Mref*a^2*epn*kpare*xm*xx*ym*cc*ss*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20)-...
            (8*3^(1-epn)*Mref*a^2*epn*kpare*xm*xx*yy*cc*ss*(-(2*cy*sx-20)/...
            Mref)^epn*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2))/(2*cy*sx-20))/...
            (9*Mref^2*xx*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)^2)-...
             Mref*pcourr*cc*((2*a*sx*(xm-xx)*(10*cy+sx+2*sy-2*cy2*sx))/...
             (3*Mref*xx*((xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^...
             (1/2))+(4*a*cx*(ym-yy)*(cy-5*sy+cy*ss))/(3*Mref*xx*((xm^2-2*...
             ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^(1/2)))+(a*csie*(2*xm^4*...
             cx*cy+12*xx^4*cx*cy-10*xm^4*cx*sy-60*xx^4*cx*sy-2*xm*ym^3*ss-...
             2*xm^3*ym*ss+4*xx*ym^3*ss+4*xx^3*ym*ss+2*xm*yy^3*ss+...
             2*xm^3*yy*ss-4*xx*yy^3*ss-4*xx^3*yy*ss+22*xm^2*xx^2*cx*cy+2*...
             xm^2*ym^2*cx*cy+10*xx^2*ym^2*cx*cy+2*xm^2*yy^2*cx*cy+...
             10*xx^2*yy^2*cx*cy-110*xm^2*xx^2*cx*sy-10*xm^2*ym^2*cx*sy-50*...
             xx^2*ym^2*cx*sy-10*xm^2*yy^2*cx*sy-50*xx^2*yy^2*cx*sy+...
             xm*ym^3*cy2*sx2+xm^3*ym*cy2*sx2-2*xx*ym^3*cy2*sx2-2*xx^3*ym*...
             cy2*sx2-xm*yy^3*cy2*sx2-xm^3*yy*cy2*sx2+2*xx*yy^3*cy2*sx2+2*...
             xx^3*yy*cy2*sx2-xm*ym^3*sx2*sy2-xm^3*ym*sx2*sy2+2*xx*ym^3*sx2*...
             sy2+2*xx^3*ym*sx2*sy2+xm*yy^3*sx2*sy2+xm^3*yy*sx2*sy2-2*xx*...
             yy^3*sx2*sy2-2*xx^3*yy*sx2*sy2-12*a*xx^5*cy*sx-24*xm*xx^3*cx*...
             cy-10*xm^3*xx*cx*cy+60*a*xx^5*ss+120*xm*xx^3*cx*sy+50*xm^3*xx*...
             cx*sy-10*xm*ym^3*cy*sx-10*xm^3*ym*cy*sx+20*xx*ym^3*cy*sx+20*...
             xx^3*ym*cy*sx+10*xm*yy^3*cy*sx+10*xm^3*yy*cy*sx-20*xx*yy^3*...
             cy*sx-20*xx^3*yy*cy*sx+4*xm*xx^2*ym*cy2*sx2-4*xm^2*xx*ym*cy2*...
             sx2-4*xm*xx^2*yy*cy2*sx2+4*xm^2*xx*yy*cy2*sx2+3*xm*ym*yy^2*...
             cy2*sx2-3*xm*ym^2*yy*cy2*sx2-6*xx*ym*yy^2*cy2*sx2+6*xx*ym^2*...
             yy*cy2*sx2-4*xm*xx^2*ym*sx2*sy2+4*xm^2*xx*ym*sx2*sy2+4*xm*...
             xx^2*yy*sx2*sy2-4*xm^2*xx*yy*sx2*sy2-3*xm*ym*yy^2*sx2*sy2+3*...
             xm*ym^2*yy*sx2*sy2+6*xx*ym*yy^2*sx2*sy2-6*xx*ym^2*yy*sx2*sy2+...
             40*a*xx^4*ym*cx*cy-40*a*xx^4*yy*cx*cy+20*a*xm*xx^4*cy*sx-2*a*...
             xm^4*xx*cy*sx-2*a*xx*ym^4*cy*sx+8*a*xx^4*ym*cx*sy-2*a*xx*yy^4*...
             cy*sx-8*a*xx^4*yy*cx*sy-6*xm*xx*ym^2*cx*cy-6*xm*xx*yy^2*cx*cy-...
             4*xm^2*ym*yy*cx*cy-20*xx^2*ym*yy*cx*cy-100*a*xm*xx^4*ss+10*a*...
             xm^4*xx*ss+10*a*xx*ym^4*ss+10*a*xx*yy^4*ss+30*xm*xx*ym^2*cx*sy-...
             40*xm*xx^2*ym*cy*sx+40*xm^2*xx*ym*cy*sx+30*xm*xx*yy^2*cx*sy+40*...
             xm*xx^2*yy*cy*sx-40*xm^2*xx*yy*cy*sx-30*xm*ym*yy^2*cy*sx+30*xm*...
             ym^2*yy*cy*sx+20*xm^2*ym*yy*cx*sy+60*xx*ym*yy^2*cy*sx-60*xx*...
             ym^2*yy*cy*sx+100*xx^2*ym*yy*cx*sy-8*xm*xx^2*ym*ss+8*xm^2*xx*...
             ym*ss+8*xm*xx^2*yy*ss-8*xm^2*xx*yy*ss-6*xm*ym*yy^2*ss+6*xm*...
             ym^2*yy*ss+12*xx*ym*yy^2*ss-12*xx*ym^2*yy*ss+20*a*xx^2*ym^3*...
             cx*cy-20*a*xx^2*yy^3*cx*cy-18*a*xm^2*xx^3*cy*sx+8*a*xm^3*xx^2*...
             cy*sx+4*a*xx^2*ym^3*cx*sy-10*a*xx^3*ym^2*cy*sx-4*a*xx^2*yy^3*...
             cx*sy-10*a*xx^3*yy^2*cy*sx+90*a*xm^2*xx^3*ss-40*a*xm^3*xx^2*...
             ss+2*xm^4*cc*ss+12*xx^4*cc*ss+50*a*xx^3*ym^2*ss+50*a*xx^3*yy^2*...
             ss+8*a*xx^5*cx2*cy*sy-16*a*xx^5*cy*sx2*sy+60*a*xm^2*xx^2*ym*...
             cx*cy-60*a*xm^2*xx^2*yy*cx*cy+60*a*xx^2*ym*yy^2*cx*cy-60*a*...
             xx^2*ym^2*yy*cx*cy+8*a*xm*xx^2*ym^2*cy*sx-4*a*xm^2*xx*ym^2*...
             cy*sx+12*a*xm^2*xx^2*ym*cx*sy+8*a*xm*xx^2*yy^2*cy*sx-4*a*xm^2*...
             xx*yy^2*cy*sx-12*a*xm^2*xx^2*yy*cx*sy-12*a*xx*ym^2*yy^2*cy*...
             sx+12*a*xx^2*ym*yy^2*cx*sy-12*a*xx^2*ym^2*yy*cx*sy-24*xm*xx^3*...
             cc*ss-10*xm^3*xx*cc*ss-40*a*xm*xx^2*ym^2*ss+20*a*xm^2*xx*ym^2*...
             ss-40*a*xm*xx^2*yy^2*ss+20*a*xm^2*xx*yy^2*ss+60*a*xx*ym^2*yy^2*...
             ss-8*a*xx^4*ym*cx*cy2*sx-16*a*xm*xx^4*cx2*cy*sy+2*a*xm^4*xx*cx2*...
             cy*sy+8*a*xx^4*yy*cx*cy2*sx+8*a*xx^4*ym*cx*sx*sy2+24*a*xm*xx^4*...
             cy*sx2*sy-2*a*xm^4*xx*cy*sx2*sy-8*a*xx^4*yy*cx*sx*sy2-4*a*xx*...
             ym^4*cy*sx2*sy-4*a*xx*yy^4*cy*sx2*sy+12*xm*xx*ym*yy*cx*cy-60*...
             xm*xx*ym*yy*cx*sy+22*xm^2*xx^2*cc*ss+2*xm^2*ym^2*cc*ss+10*xx^2*...
             ym^2*cc*ss+2*xm^2*yy^2*cc*ss+10*xx^2*yy^2*cc*ss-4*a*xx^2*ym^3*...
             cx*cy2*sx+16*a*xm^2*xx^3*cx2*cy*sy-8*a*xm^3*xx^2*cx2*cy*sy+4*...
             a*xx^2*yy^3*cx*cy2*sx+4*a*xx^3*ym^2*cx2*cy*sy+4*a*xx^3*yy^2*cx2*...
             cy*sy-20*a*xm*xx*ym^3*cx*cy-80*a*xm*xx^3*ym*cx*cy-20*a*xm^3*xx*...
             ym*cx*cy+20*a*xm*xx*yy^3*cx*cy+80*a*xm*xx^3*yy*cx*cy+20*a*xm^3*...
             xx*yy*cx*cy+4*a*xx^2*ym^3*cx*sx*sy2-20*a*xm^2*xx^3*cy*sx2*sy+8*...
             a*xm^3*xx^2*cy*sx2*sy-4*a*xx^2*yy^3*cx*sx*sy2-16*a*xx^3*ym^2*cy*...
             sx2*sy-16*a*xx^3*yy^2*cy*sx2*sy-4*a*xm*xx*ym^3*cx*sy-16*a*xm*...
             xx^3*ym*cx*sy-4*a*xm^3*xx*ym*cx*sy+4*a*xm*xx*yy^3*cx*sy+16*a*...
             xm*xx^3*yy*cx*sy+4*a*xm^3*xx*yy*cx*sy+8*a*xx*ym*yy^3*cy*sx+8*...
             a*xx*ym^3*yy*cy*sx+20*a*xx^3*ym*yy*cy*sx-40*a*xx*ym*yy^3*ss-40*...
             a*xx*ym^3*yy*ss-100*a*xx^3*ym*yy*ss-12*a*xm^2*xx^2*ym*cx*cy2*...
             sx+12*a*xm^2*xx^2*yy*cx*cy2*sx-4*a*xm*xx^2*ym^2*cx2*cy*sy+2*...
             a*xm^2*xx*ym^2*cx2*cy*sy-12*a*xx^2*ym*yy^2*cx*cy2*sx+12*a*xx^2*...
             ym^2*yy*cx*cy2*sx-4*a*xm*xx^2*yy^2*cx2*cy*sy+2*a*xm^2*xx*yy^2*...
             cx2*cy*sy-60*a*xm*xx*ym*yy^2*cx*cy+60*a*xm*xx*ym^2*yy*cx*cy+...
             12*a*xm^2*xx^2*ym*cx*sx*sy2-12*a*xm^2*xx^2*yy*cx*sx*sy2+12*a*...
             xm*xx^2*ym^2*cy*sx2*sy-6*a*xm^2*xx*ym^2*cy*sx2*sy+12*a*xx^2*...
             ym*yy^2*cx*sx*sy2-12*a*xx^2*ym^2*yy*cx*sx*sy2+12*a*xm*xx^2*...
             yy^2*cy*sx2*sy-6*a*xm^2*xx*yy^2*cy*sx2*sy-24*a*xx*ym^2*yy^2*...
             cy*sx2*sy-12*a*xm*xx*ym*yy^2*cx*sy+12*a*xm*xx*ym^2*yy*cx*sy-...
             16*a*xm*xx^2*ym*yy*cy*sx+8*a*xm^2*xx*ym*yy*cy*sx+80*a*xm*xx^2*...
             ym*yy*ss-40*a*xm^2*xx*ym*yy*ss-6*xm*xx*ym^2*cc*ss-6*xm*xx*yy^2*...
             cc*ss-4*xm^2*ym*yy*cc*ss-20*xx^2*ym*yy*cc*ss+4*a*xm*xx*ym^3*...
             cx*cy2*sx+16*a*xm*xx^3*ym*cx*cy2*sx+4*a*xm^3*xx*ym*cx*cy2*sx-...
             4*a*xm*xx*yy^3*cx*cy2*sx-16*a*xm*xx^3*yy*cx*cy2*sx-4*a*xm^3*...
             xx*yy*cx*cy2*sx-8*a*xx^3*ym*yy*cx2*cy*sy-4*a*xm*xx*ym^3*cx*sx*...
             sy2-16*a*xm*xx^3*ym*cx*sx*sy2-4*a*xm^3*xx*ym*cx*sx*sy2+4*a*xm*...
             xx*yy^3*cx*sx*sy2+16*a*xm*xx^3*yy*cx*sx*sy2+4*a*xm^3*xx*yy*cx*...
             sx*sy2+16*a*xx*ym*yy^3*cy*sx2*sy+16*a*xx*ym^3*yy*cy*sx2*sy+32*...
             a*xx^3*ym*yy*cy*sx2*sy+12*a*xm*xx*ym*yy^2*cx*cy2*sx-12*a*xm*xx*...
             ym^2*yy*cx*cy2*sx+8*a*xm*xx^2*ym*yy*cx2*cy*sy-4*a*xm^2*xx*ym*...
             yy*cx2*cy*sy-12*a*xm*xx*ym*yy^2*cx*sx*sy2+12*a*xm*xx*ym^2*yy*...
             cx*sx*sy2-24*a*xm*xx^2*ym*yy*cy*sx2*sy+12*a*xm^2*xx*ym*yy*cy*...
             sx2*sy+12*xm*xx*ym*yy*cc*ss))/(xx*(xm^2-2*ym*yy-2*xm*xx+2*xx^2+...
             ym^2+yy^2)^2)-(3^(1/2)*(ss+2)^2*(2*cx*sy-cx2*cy2+2*cy*sx+20))/...
             (2*tie*(cy*sx-10)*(-(2*cy*sx-20)/Mref)^(1/2))+(5*a*cx*sy*(cy*sx-10)*...
             (ss+2)*(xm-xx))/(3*xx*((xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/...
             xx^2)^(1/2))+(5*a*cc*sx*(xm-xx)*(10*cy+sx+2*sy-2*cy2*sx))/(3*xx*...
             ((xm^2-2*ym*yy-2*xm*xx+2*xx^2+ym^2+yy^2)/xx^2)^(1/2))-(5*cc*(cy*...
             sx-10)*(2*ym-2*yy)*(ss+2)*(xm-xx))/(6*xx^3*((xm^2-2*ym*yy-2*xm*...
             xx+2*xx^2+ym^2+yy^2)/xx^2)^(3/2)))
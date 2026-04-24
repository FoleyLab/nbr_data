                    tv = float(theta_match.group(1))
                    if tv > 100:
                        tv_flip = 180 - tv
                        tv = tv_flip
                    theta_vals.append(tv)

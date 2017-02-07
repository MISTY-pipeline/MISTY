import yt

## getting useful hidden functions from AbsorptionSpectrum in yt


##field_data = ray.all_data()
## use_peculiar_velocity is a boolean for yes/no
## comments added by JT 


def apply_observing_redshift(self, field_data, use_peculiar_velocity,
                                 observing_redshift):
        """
        Change the redshifts of individual absorbers to account for the
        redshift at which the observer sits.

        The intermediate redshift that is seen by an observer
        at a redshift other than z=0 is z12, where z1 is the
        observing redshift and z2 is the emitted photon's redshift
        Hogg (2000) eq. 13:

        1 + z12 = (1 + z2) / (1 + z1)
        """
        if observing_redshift == 0.:
            # This is already assumed in the generation of the LightRay
            redshift = field_data['redshift']
            if use_peculiar_velocity:
                redshift_eff = field_data['redshift_eff']
        else:
            # The intermediate redshift that is seen by an observer
            # at a redshift other than z=0 is z12, where z1 is the
            # observing redshift and z2 is the emitted photon's redshift
            # Hogg (2000) eq. 13:
            # 1 + z12 = (1 + z2) / (1 + z1)
            redshift = ((1 + field_data['redshift']) / \
                        (1 + observing_redshift)) - 1.
            # Combining cosmological redshift and doppler redshift
            # into an effective redshift is found in Peacock's
            # Cosmological Physics eqn 3.75:
            # 1 + z_eff = (1 + z_cosmo) * (1 + z_doppler)
            if use_peculiar_velocity:
                redshift_eff = ((1 + redshift) * \
                                (1 + field_data['redshift_dopp'])) - 1.

        if not use_peculiar_velocity:
            redshift_eff = redshift

        return redshift, redshift_eff



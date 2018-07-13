import os
from conan.packager import ConanMultiPackager


if __name__ == "__main__":
    # ConanPackageTools
    # See: https://github.com/conan-io/conan-package-tools/blob/develop/README.md
    builder = ConanMultiPackager()
    builder.add_common_builds()
    updated_builds = []
    # Add environment variables to build definitions
    xms_version = os.environ.get('XMS_VERSION', None)

    for settings, options, env_vars, build_requires, reference in builder.items:
        env_vars.update({
            'XMS_VERSION': xms_version,
            'VERBOSE': 1
        })

        # Require c++11 compatibility
        if settings['compiler'] == 'gcc':
            settings.update({
                'compiler.libcxx': 'libstdc++11'
            })

        # Add additional configurations for xms options
        if settings['compiler'] == 'Visual Studio' and 'MD' in settings['compiler.runtime']:
            xms_options = dict(options)
            xms_options.update({'xmsmesh:xms': True})
            updated_builds.append([settings, xms_options, env_vars, build_requires])

        updated_builds.append([settings, options, env_vars, build_requires])

    builder.builds = updated_builds

    builder.run()

#include "cli.h"

#include <exception>
#include <iostream>

#include <CLI/CLI.hpp>

#include "app/subcommands.h"

namespace
{
template <typename Function>
int run_command(Function run)
{
    try
    {
        run();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
}  // namespace

namespace vcfbox
{
int run_cli(int argc, char** argv)
{
    CLI::App app{"A simple cli tool contains some operation on vcf"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(1);

    CombineOptions combine_options;
    ParentageOptions parentage_options;

    auto* combine_command = add_combine_command(app, combine_options);
    auto* parentage_command = add_parentage_command(app, parentage_options);

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    if (*combine_command)
    {
        return run_command([&] { run_combine_command(combine_options); });
    }
    if (*parentage_command)
    {
        return run_command([&] { run_parentage_command(parentage_options); });
    }

    return 0;
}
}  // namespace vcfbox

{% set fullname = fullname | replace(module ~ ".", "")%}
{% if objtype == "method" %}
    {% set fullname = fullname ~ "()" %}
{% endif %}

{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

Lives in: {{ module }}

.. auto{{ objtype }}:: {{ objname }}